/******* libnumeric.c *******//*
Copyright (C) 2014-2015 Ivan Markin
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

#include <math.h>
#include <complex.h>
#include <malloc.h>
#include "../libmesh/libmesh.h"
#include <fftw3.h>
#include "libnumeric.h"


double mass_of(mesh * space, double * rho) {
	/* Initial mass */
	double mass = 0.0;
	/* Integrating mass */
	for (point j=0; j< space->points; j++)
		mass += 4*PI*pow(space->map[j],2)*rho[j] * space->avg_res;
	return mass;
}

int solve_tridiagonal_sweep_inplace( Complex *a, Complex *b, Complex *c, Complex *d, Complex *x, equation n_eqs) {
	/*= Solving =*/
	printf("libnumeric: [i] Sovling started\n");
	/* Initiation */
	c[0]/=b[0];
	d[0]/=b[0];
	
	/* Forward speed */
	printf("libnumeric: [i] Forward sweep\n");
	for (equation i=1; i<n_eqs; i++) {
		c[i]=c[i]/(b[i]-c[i-1]*a[i]);
		d[i]=(d[i]-(d[i-1]*a[i]))/(b[i]-(c[i-1]*a[i]));
	}
	/* Backward sweep */
	printf("libnumeric: [i] Backward sweep\n");
	x[n_eqs-1] = d[n_eqs-1];
	equation i=n_eqs-1;
	do {	i--;
		x[i]=d[i]-c[i]*x[i+1];
	} while(i!=0);

	return 0;
}

int solve_poisson_sweep(mesh * space, Complex * U, double * rho){
	printf("libnumeric: [!] Hello from C world!\n");
	/* Session constants */
	double res_sq = space->avg_res*space->avg_res; 
	double M = 1.; // mass_of(space, rho);
	equation n_eqs = space->points-1; // Number of equations
	/* Matrix creation */		
	Complex *a,*b,*c, *d;
	a = (Complex *) malloc( n_eqs*sizeof(Complex));
	b = (Complex *) malloc( n_eqs*sizeof(Complex));
	c = (Complex *) malloc( n_eqs*sizeof(Complex));

	d = (Complex *) malloc( n_eqs*sizeof(Complex));
	/* Boundary conditions */
	/*- left -*/
	a[0] = 0.; // Always 0
	b[0] = -2.;
	c[0] = 2.;
	d[0] = res_sq * rho[0];
	/*- right -*/
	U[space->LAST] = -M/space->map[space->LAST];
	
	a[n_eqs-1] = 1 - space->avg_res/space->map[n_eqs-1];
	b[n_eqs-1] = -2.;
	c[n_eqs-1] = 0.;
	d[n_eqs-1] = res_sq * rho[n_eqs-1]  - (1 + space->avg_res/(space->map[n_eqs-1])) * U[space->LAST];
	/* Filling matrix */
	printf("libnumeric: [i] Filling matrix...\n");
	for (equation j=1; j < n_eqs-1; j++) {
		a[j] = 1 - space->avg_res/space->map[j];
		b[j] = -2.;
		c[j] = 1 + space->avg_res/space->map[j];

		d[j] = res_sq * rho[j];
	}	
	
	solve_tridiagonal_sweep_inplace(a, b, c, d, U, n_eqs);	
	U[0] = 3.*U[1] - 3.*U[2] + U[3];
	/* Freeing matrix */
	free(a);
	free(b);
	free(c);
	
	free(d);	
	/* Bye */
	printf("libnumeric: [+] Solving done!\n");
	return 1;
}
int solve_poisson_sweep_convar(mesh * space, Complex * V, double * rho){
	printf("libnumeric: [!] Hello from C world!\n");
	/* Session constants */
	double res_sq = space->avg_res*space->avg_res; 
	double M = 1.; // mass_of(space, rho);
	equation n_eqs = space->points-2; //sure??// Number of equations
	
	V[space->LAST] = -M;
	/* Matrix creation */		
	Complex *a,*b,*c, *d;
	a = (Complex *) malloc( n_eqs *sizeof(Complex));
	b = (Complex *) malloc( n_eqs *sizeof(Complex));
	c = (Complex *) malloc( n_eqs *sizeof(Complex));
	d = (Complex *) malloc( n_eqs *sizeof(Complex));
	/* Boundary conditions */
	/*- left -*/
	a[0] = 0.;
	b[0] = -2.;
	c[0] = 1.;
	d[0] = res_sq * rho[1];
	/*- right -*/
	a[n_eqs-1] = 1.;
	b[n_eqs-1] = -2.;
	c[n_eqs-1] = 0.;
	d[n_eqs-1] = res_sq * rho[space->points-2] - V[space->LAST];
	/* Filling matrix */
	printf("libnumeric: [i] Filling matrix...\n");
	for (equation j=1; j < n_eqs-1; j++) {
		a[j] = 1.;
		b[j] = -2.;
		c[j] = 1.;
		d[j] = res_sq * rho[j+1];
	}	
	
	solve_tridiagonal_sweep_inplace(a, b, c, d, V, n_eqs);	
	

	V[0] = 3.*V[1] - 3.*V[2] + V[3];
	
	/* Freeing matrix */
	free(a);
	free(b);
	free(c);
	free(d);	
	/* Bye */
	printf("libnumeric: [+] Solving done!\n");
	return 1;
}

/*
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
				psi_new[j] = psi[j] - M*I* (  (psi_old[j+1] + psi_old[j-1 + j?0:2 ] - 2.*psi_old[j] ) + ( psi[j+1] + psi[j-1 + j?0:2] -2.*psi[j] ) + j?(1/(2.*space->res*space->map[j]) *  ( (psi_old[j+1]-psi_old[j-1]) + (psi[j+1]-psi[j-1])  ) ):0  ) -  I * T2 * V[j] * ( psi_old[j] + psi[j] );
			for(int j=1; j< space->points-1; j++) 
				psi_old[j]=psi_new[j]; 
			
		} // Now we have solution for n'th psi
	}

	return 1;
}
*/
