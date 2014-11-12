/******* libnumeric.c *******//*
Copyright (C) 2014 Ivan Markin
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>. */

//#include "libnumeric.h"
#include <math.h>
#include <complex.h>
#include <malloc.h>
#include "../libmesh/libmesh.h"
#include <fftw3.h>

//#pragma STDC CS_LIMITED_RANGE on

#define PI 3.141592653589793238462643
#define dPI 6.2831853071795864769252867 // 2*PI
#define A 0.398942280401432677939946 // 1/sqrt(2*PI)

#define MAX_CRANK_NICKOLSON_ITER 3

typedef double complex Complex; 
typedef int points;

//typedef double dim;
typedef unsigned long int dot;


double mass_of(mesh * space, double * rho) {
	double mass = 0.0;
	
	for (dot j=0; j< space->points; j++) {
		mass += rho[j] * space->res;
		//printf("d_mass: %lf\n", rho[j]);
	}
	return mass;
}


int solve_poisson_sweep(mesh * space, Complex * U, double * rho){
	//for (int j=0; j< space->points; j++)
	//	printf("rhO: %lf\n", rho[j]);

	double dx2 = space->res*space->res; //Const
			
	Complex *a,*b,*c, *d;
	a = (Complex *) malloc(space->points*sizeof(Complex));
	b = (Complex *) malloc(space->points*sizeof(Complex));
	c = (Complex *) malloc(space->points*sizeof(Complex));

	d = (Complex *) malloc(space->points*sizeof(Complex));
	
	// Boundary conditions
	double M = mass_of(space, rho);
	//printf("Mass: %lf\n", M);
	a[0] = 0.; b[0] = -2.; c[0] = 2.;	d[0] = dx2 * rho[0];
	a[space->points - 1] = 1 - space->res/(2.*space->map[space->points-1]);
	b[space->points - 1] = -2.;
	c[space->points - 1] = 0.;
	d[space->points - 1] = dx2 * rho[space->points-1]  - (1 + space->res/(2.*space->map[space->points-1])*M/(space->map[space->points-1] + space->res));

	// Filling matrix
		for (dot j=1; j < space->points - 1; j++) {
			
			a[j] = (1 - space->res/(2.*space->map[j]));
			b[j] = 2.;
			c[j] = (1 + space->res/(2.*space->map[j]));

			d[j]=dx2 * rho[j];
		}	

		// Starting sovle
		c[0]=c[0]/b[0];
		d[0]=d[0]/b[0];

		//Forward sweep
		for (dot i=1; i<space->points; i++) {
			c[i]=c[i]/(b[i]-c[i-1]*a[i]);
		
			d[i]=(d[i]-(d[i-1]*a[i]))/(b[i]-(c[i-1]*a[i]));
		}
		//Backward sweep
		U[space->points-1]=d[space->points-1];
		for (dot i=space->points-1; i!=0; --i) {
			U[i]=d[i]-c[i]*U[i+1];
		}

	//for (int j=0; j< space->points; j++)
	//	printf("rho after poisson: %lf\n", rho[j]);


	free(a);
	free(b);
	free(c);
	
	free(d);	

	return 1;
}


int solve_spherically_symmetric(mesh * space, mesh * time, Complex * psi) {
	
	double M = time->res/(4.*space->res*space->res);
	double T2 = time->res/2;
	
	Complex * psi_new; //Wave function current approximation
	Complex * psi_old = psi;
	psi_new = (Complex *) malloc (space->points * sizeof(Complex));
	
	//double * rho = (double *) malloc(space->points * sizeof(double));

	Complex * V = (Complex *) malloc (space->points * sizeof(Complex)); // Potential
	for (int j=0; j < space->points; j++) {
		V[j] =0.0;
	}

	for (int n=0; n< time->points; n++) {
	 	//printf("\ntime %d: ", n);
		// Copying old to new as 0'th approx.
		for (int i=0; i<space->points; i++)
			psi_new[i]=psi_old[i];
		//Starting iterations	
		for(int it=0; it < MAX_CRANK_NICKOLSON_ITER; it++) {
			
			//for (int j=0; j< space->points; j++)
			//	rho[j] = fabs(psi_new[j]) - 1.; //Should be average
			
			// Getting solution for V potential
			//solve_poisson_sweep(space, V, rho); // for U
			//for( int j=0; j< space->points; 
			// Now we need to find  solution for psi at next iteration
			for(int j=0; j< (space->points - 1); j++) 
				psi_new[j] = psi[j] - M*I* (  (psi_old[j+1] + psi_old[j-1 + j?0:2 ] - 2.*psi_old[j] ) + ( psi[j+1] + psi[j-1 + j?0:2] -2.*psi[j] ) + j?(1/(4.*space->res*space->map[j]) *  ( (psi_old[j+1]-psi_old[j-1]) + (psi[j+1]-psi[j-1])  ) ):0  ) -  I * T2 * V[j] * ( psi_old[j] + psi[j] );
			for(int j=1; j< space->points-1; j++) 
				psi_old[j]=psi_new[j]; 
			
		} // Now we have solution for n'th psi
	}

	return 1;
}
