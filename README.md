libnumeric
==========

C library with some useful numerical algos.

libnumeric.c -  Solving Equations Numerically

In this library:

1 .Sweep method for tridiagonal equation

2. FFT method for Poisson equation

3. Iterative Crank-Nickolson

Compile as:

$ gcc -std=c99 -o libnumeric libnumeric.c -lm -lmesh -lconfig -lfftw3

This edition by default uses libnumeric.cfg configuration file in working directory.
You can specify configuration file in first argument like: $ ./libnumeric path/to/some/file.cfg

Use graph.gnuplot script to plot some results. Run this from working directory:

$ plot/graph.gnuplot

Results will be stored in "results" directory in csv format.
