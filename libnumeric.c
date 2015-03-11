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

int _solve_tridiagonal_sweep_inplace(trimatrix_t * T, Complex *x, equation n_eqs) {
	/* Initiation */
	T[0].c/=T[0].b;
	T[0].d/=T[0].b;
	
	/* Forward speed */
	for (equation i=1; i<n_eqs; i++) {
		T[i].c = T[i].c/(T[i].b-T[i-1].c*T[i].a);
		T[i].d = (T[i].d-(T[i-1].d*T[i].a))/(T[i].b-(T[i-1].c*T[i].a));
	}
	/* Backward sweep */
	x[n_eqs-1] = T[n_eqs-1].d;
	equation i=n_eqs-1;
	do {	i--;
		x[i] = T[i].d-T[i].c*x[i+1];
	} while(i!=0);

	return 0;
}

int solve_poisson_sweep(mesh * space, Complex * U, double * rho){
	/* Session constants */
	double res_sq = space->avg_res*space->avg_res; 
	double M = 1.; // mass_of(space, rho);
	equation n_eqs = space->points-1; //sure??// Number of equations
	U[space->LAST] = -M/space->map[space->LAST];
	
	/* Matrix creation */		
	trimatrix_t *T = (trimatrix_t *) malloc( n_eqs * sizeof(trimatrix_t));
	/* Boundary conditions */
	/*- left -*/
	T[0] = (trimatrix_t) { 0., -2., 2., res_sq * rho[0] };
	/*- right -*/
	T[n_eqs-1] = (trimatrix_t) { 1.- space->avg_res/space->map[n_eqs-1], -2., 0., res_sq * rho[n_eqs-1] - (1 + space->avg_res/(space->map[n_eqs-1])) * U[space->LAST] };
	for (equation j=1; j < n_eqs-1; j++) 
		T[j] = (trimatrix_t){ 1. - space->avg_res/space->map[j], -2., 1. + space->avg_res/space->map[j], res_sq * rho[j] };
	
	_solve_tridiagonal_sweep_inplace(T, U, n_eqs);	

	U[0] = 3.*U[1] - 3.*U[2] + U[3];
	/* Freeing matrix */
	free(T);
	return 1;
}
int solve_poisson_sweep_convars(mesh * space, Complex * V, double * rho){
	/* Session constants */
	double res_sq = space->avg_res*space->avg_res; 
	double M = 1.; // mass_of(space, rho);
	equation n_eqs = space->points-2; //sure??// Number of equations
	
	V[space->LAST] = -M;
	/* Matrix creation */		
	trimatrix_t *T = (trimatrix_t *) malloc( n_eqs * sizeof(trimatrix_t));
	/* Boundary conditions */
	/*- left -*/
	T[0] = (trimatrix_t) { 0., -2., 1., res_sq * rho[1] };
	/*- right -*/
	T[n_eqs-1] = (trimatrix_t) { 1., -2., 0., res_sq * rho[space->points-2] - V[space->LAST] };
	/* Filling matrix */
	for (equation j=1; j < n_eqs-1; j++) 
		T[j] = (trimatrix_t){ 1., -2., 1., res_sq * rho[j+1] };
	
	_solve_tridiagonal_sweep_inplace(T, V+1, n_eqs);	

	V[0] = 0.;
	
	/* Freeing matrix */
	free(T);
	return 1;
}

int solve_shroedinger_sweep_convars(mesh * space, mesh * time, Complex * phi, Complex * V, Complex psi_0 ){
	/* Session constants */
	double s_res_sq = space->avg_res*space->avg_res; 
	double s2_t_res = 4.*s_res_sq/time->res;
	equation n_eqs = time->points-2;
	
	/* Matrix creation */		
	trimatrix_t *T = (trimatrix_t *) malloc( n_eqs * sizeof(trimatrix_t));
	/* Filling matrix */
	for (int j=1; j <= n_eqs; j++) 
		T[j-1] = (trimatrix_t){ 1., 2.(1.- s_res_sq*V[j]/space->map[j]), 1. + I*s2_t_res,
					I*s2_t_res*phi[j] - (phi[j-1] -2.*phi[j] + phi[j+1]) + 2.*s_res_sq * V[j]/space.map[j] * phi[j] + 4.*s_res_sq*V[j]*psi_0[j] };
	/* Boundary conditions */
	/*- left -*/
	T[0].d -= T[0].a * phi[0];
	T[0].a = 0.;
	/*- right -*/
	T[n_eqs-1].d -= T[n_eqs-1].c * phi[time->points];
	T[n_eqs-1].c = 0.;
	
	_solve_tridiagonal_sweep_inplace(T, phi+1, n_eqs);	
	
	/* Freeing matrix */
	free(T);
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
