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

Optional  
`libconfig` - library for pretty looking configuration files for `solver.c`.    
Building
--------

Complie as:  

`$ gcc -std=c99 -o libnumeric libnumeric.c -lm -lmesh -lconfig -lfftw3`  

This edition by default uses `libnumeric.cfg` configuration file in working directory.  
You can specify configuration file in first argument like:  
`$ ./libnumeric path/to/some/file.cfg`  

Plotting results
--------

Use graph.gnuplot script to plot some results. Run this from working directory:  

`$ plot/graph.gnuplot`  

Results will be stored in `results` directory in csv format.
