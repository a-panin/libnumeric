/*	libnumeric.c -  Solving Equations Numerically		*/
/*	In this library: 					*/
/*	1 .Sweep method for tridiagonal equation		*/
/*	2. FFT method for Poisson equation			*/
/*	3. Iterative Crank-Nickolson				*/

/*	Compile as:								*/
/*	$ gcc -std=c99 -o libnumeric libnumeric.c -lm -lmesh -lconfig -lfftw3	*/

/*	This edition by default uses libnumeric.cfg configuration file in working directory	*/
/*	You can specify configuration file in first argument like: $ ./libnumeric path/to/some/file.cfg	*/

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <malloc.h>
#include <mesh.h>
#include <libconfig.h>
#include <fftw3.h>
#pragma STDC CS_LIMITED_RANGE on

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
	double temp = (x-x0)/(2.*sigma); 
	double ssigma = sigma*sigma;
	double pp0 = p0*p0;
	Complex c = A*sigma;
	c = cpow(c,0.5);
	c = c * (exp(-ssigma*pp0));
	
	Complex c2 = ssigma+I*t/2.;
	c2 = cpow(c2,-0.5);
	
	c = c*c2;
	
	c2 = x-2.*ssigma*p0*I;
	c2 = (-1.)*c2*c2;

	Complex c3 = 4.*ssigma + 2.*t*I;
		
	return c*cexp(c2/c3);
}



Complex packet_framework(double x, int n, double * p) {
	if (n >=0 && n<=4) p[n]=x;	
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

// End of shrodinger sweep

// Start of poisson
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
	for (int k=1;k< space->points; k++) {
		C[k]=dS/(cos(cC*k)-1);
	}
	
	fftw_plan p;
	p = fftw_plan_dft_1d(space->points,rho,tmp,FFTW_FORWARD,FFTW_ESTIMATE); //Use FFTW_MEASURE in Release, because it should run faster.
 
	fftw_execute(p);
	
	for (int i=0; i<space->points; i++) {
		tmp[i]=tmp[i]*C[i]/space->points;
	}


	p = fftw_plan_dft_1d(space->points,tmp,u,FFTW_BACKWARD, FFTW_ESTIMATE); 

	fftw_execute(p);

	fftw_destroy_plan(p);
	free(tmp);
	free(C);

	return 1;
}
/// Poisson end


/// Iterative Crank-Nickolson start

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


// ICN end


int main(int argc, char *argv[]){

	Complex * psi; // Wave function snapshot
	mesh space = mesh_default;
	mesh time = mesh_default;
	double p[5] = {0., 0., 0., 0., 0.}; // Packet initial parameters
	
	/// Config loading 
	config_t config;
	config_init(&config);
	if (config_read_file(&config, (argc==1) ? "libnumeric.cfg" : argv[1] )==CONFIG_FALSE) {
		printf("Configuration problem: %s\n",config_error_text(&config));		
		return -1;
	}
	
	if (	
	config_lookup_float(&config,"mesh.space.left",&space.left) *
	config_lookup_float(&config,"mesh.space.right",&space.right) *	
	config_lookup_int(&config,"mesh.space.points",&space.points) *	
	
	config_lookup_float(&config,"mesh.time.left",&time.left) *	
	config_lookup_float(&config,"mesh.time.right",&time.right) *	
	config_lookup_int(&config,"mesh.time.points",&time.points) *
	
	config_lookup_float(&config,"init.gauss_packet.mean",&p[2]) *	
	config_lookup_float(&config,"init.gauss_packet.sigma",&p[3]) *	
	config_lookup_float(&config,"init.gauss_packet.p0",&p[4]) 
	== 0
	) {
		printf("Config parameters problem! Try heck points in floats!\n");
		return -1;
	}
	config_destroy(&config);
	// Config loaded


	recalc(&space);
	recalc(&time);


	psi = (Complex *) malloc(space.points*sizeof(Complex)); //Include into recalc? 
	mesh_flush(psi,&space,&packet_framework, 0, p);
	/* Choose solving method */

	//solve_shrodinger_sweep(&space, &time, psi);
	//solve_poisson_fft(&space, psi, psi);	
	solve_ICN(&space, &time, psi);
	
	/* * * * */

	printf("Done, writing result...\n");
	write_power_to_csv("results/psi.csv",psi,&space);
	
	// Rename to destroy!
	destruct(&space);
	destruct(&time);
	return 0;
}

