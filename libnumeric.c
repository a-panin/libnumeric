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

#include <math.h>
#include <complex.h>
#include <malloc.h>
#include <libmesh.h>
#include <fftw3.h>

//#pragma STDC CS_LIMITED_RANGE on

#define PI 3.141592653589793238462643
#define dPI 6.2831853071795864769252867 // 2*PI
#define A 0.398942280401432677939946 // 1/sqrt(2*PI)

#define CRANK_NICKOLSON_ITER 3

typedef double complex Complex; 
typedef int points;

//typedef double dim;
typedef unsigned long int dot;

void operation_on_array (Complex * child, Complex * parent, int size, Complex (*f) (Complex)) {
	for (int i=0; i < size; i++) {
		child[i] = (*f) ( parent[i]); 	
	}
}


//inline ?
Complex norm(Complex z) {
	return pow(cabs(z),2.);
}

Complex packet (double x, double t, double x0, double sigma, double p0 ) {
	return sqrt( (sigma/sqrt(2 * PI ) ) /  ( pow(sigma,2) + (I*t/2.) ) ) * exp(-pow(sigma*p0,2)) * exp( -0.25*pow(x-2*I*pow(sigma*p0,2),2)/(pow(sigma,2) +(I*t/2)) );  
}



Complex packet_framework(double x, int n, double * p) {
	if (n >=0 && n<=4)
		p[n]=x;	
	return packet(p[0],p[1],p[2],p[3],p[4]);
}




int solve_shrodinger_sweep(mesh * space, mesh * time, Complex * psi){
	double M = 2.*(space->res*space->res)/time->res; //Const
			
	Complex *a,*b,*c, *d;
	a = (Complex *) malloc(space->points*sizeof(Complex));
	b = (Complex *) malloc(space->points*sizeof(Complex));
	c = (Complex *) malloc(space->points*sizeof(Complex));

	d = (Complex *) malloc(space->points*sizeof(Complex));


	Complex act	= 1.;
	Complex bt	= -2.+M*I;
	Complex D	= M*I;

	for (dot j=0; j < time->points ; j++) {
		for (dot i=0; i < space->points ; i++) {
			
			a[i]=act;
			b[i]=bt;
			c[i]=act;

			d[i]=D*psi[i];
		}	
		// Boundary conditions should initial boundaries of psi :
		a[0]=a[0]*(0.+0.*I); /// Boundary conds !!! psi[0]
		c[space->points-1]=c[space->points-1]*(0.+0.*I); /// Boundary conds !!! psi[space->points-1]


		c[0]=c[0]/b[0];
		d[0]=d[0]/b[0];

		//Forward sweep
		for (dot i=1; i<space->points; i++) {
			c[i]=c[i]/(b[i]-c[i-1]*a[i]);
		
			d[i]=(d[i]-(d[i-1]*a[i]))/(b[i]-(c[i-1]*a[i]));
		}
		//Backward sweep
		psi[space->points-1]=d[space->points-1];
		for (dot i=space->points-1; i!=0; --i) {
			psi[i]=d[i]-c[i]*psi[i+1];
		}
	}

	free(a);
	free(b);
	free(c);
	
	free(d);	

	return 1;
}

int solve_poisson_fft(mesh * space, Complex * rho, Complex * u){
	
	/*
	Alternatively, if you have a C compiler (such as gcc) that supports the C99 revision of the
	ANSI C standard, you can use C’s new native complex type (which is binary-compatible
	with the typedef above). In particular, if you #include <complex.h> before <fftw3.h>,
	then fftw_complex is defined to be the native complex type and you can manipulate it
	with ordinary arithmetic (e.g. x = y * (3+4*I), where x and y are fftw_complex and I is
	the standard symbol for the imaginary unit);
	*/
	
	
	double dS=(space->res*space->res)/2.; 
	double cC = dPI/space->points; //cosine coeff
	double *C = (double *) malloc(sizeof(double) * space->points); ///Coefficient array

	Complex * tmp = (Complex *) malloc ( space->points * sizeof(Complex)); 
	//Filling
	C[0]=1.; // !? Need for figuring out
	for (int k=1;k< space->points; k++)
		C[k]=dS/(cos(cC*k)-1);
	
	fftw_plan p;
	p = fftw_plan_dft_1d(space->points,rho,tmp,FFTW_FORWARD,FFTW_ESTIMATE); //Use FFTW_MEASURE in Release, because it should run faster.
 
	fftw_execute(p);
	
	for (int i=0; i<space->points; i++)
		tmp[i]=tmp[i]*C[i]/space->points;

	p = fftw_plan_dft_1d(space->points,tmp,u,FFTW_BACKWARD, FFTW_ESTIMATE); 

	fftw_execute(p);

	fftw_destroy_plan(p);
	free(tmp);
	free(C);

	return 1;
}

int solve_ICN ( mesh * space, mesh * time, Complex * psi) {


	double M = time->res/(4.*space->res*space->res);
	double T2 = time->res/2;
	
	Complex * psi_apx; //Wave function current approximation
	psi_apx = (Complex *) malloc (space->points * sizeof(Complex));
	
	Complex * rho = (Complex *) malloc (space->points * sizeof(Complex)); // Array for rho in potential equation (Poisson) 
	Complex * V = (Complex *) malloc (space->points * sizeof(Complex)); // Potential
		

	for (int n=0; n< time->points; n++) {
	 	for (int i=0; i<space->points; i++) {
			psi_apx[i]=psi[i];
		}
	
		for(int it=0; it < CRANK_NICKOLSON_ITER; it++) {
			
			operation_on_array(rho, psi_apx, space->points, &norm);
			solve_poisson_fft(space, rho, V); 

			// Now we need to find  solution for next psi
			for(int j=1; j< (space->points - 1); j++){ // without boundary 
				psi_apx[j]=psi[j]+M*I*(psi_apx[j+1]+psi_apx[j-1]-2.*psi_apx[j]+psi[j+1]+psi[j-1]-2.*psi[j] ) - V[j] *T2*I*( psi[j]+psi_apx[j] );
			}

		} // Now we have solution for n'th psi
	
		for(int j=1; j< space->points-1; j++) {
				 psi[j]=psi_apx[j]; 
		};
	}

	return 1;
}

