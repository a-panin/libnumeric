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
1. `math` - standard C math library.  
2. `libmesh` - essential library for mesh. <https://github.com/mark-in/libmesh/>  
3. `FFTW3` - for Poisson equation solving. You could download it at <http://fftw.org/> or instal from repos of your distro.  

Building
--------
(`automake` will be available soon)

To build as shared library:

   `$ gcc -std=c99 -lm -lmesh -lfftw3 -c -Wall -Werror -fpic libnumeric.c`  
   `$ gcc -shared -o numeric.so libnumeric.o`  

Usage example:  
   `$ gcc -L/path/to/libnumeric -Wall -o out in.c -lnumeric`  

Then move `numeric.so` to `lib` folder and `libnumeric.h` to `include` folder.  
If you want you may use direct include into your project.  
