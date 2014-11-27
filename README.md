libnumeric
==========

C library with some useful numerical algorithms.  

libnumeric.c -  Solving Equations Numerically  

In this library:  

   1. Sweep method for tridiagonal equation

   2. FFT method for Poisson equation

   3. Iterative Crank-Nickolson

Requirements
---------
1. `libmesh` - essential library for mesh. <https://github.com/mark-in/libmesh/>  
2. `FFTW3` - for Poisson equation solving. You could download it at <http://fftw.org/> or instal from repos of your distro.  

Building
--------
To build as shared library:

   `$ make`  
For clean-up: 
   `$ make clean`  

Now you have shared object.
