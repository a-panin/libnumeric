CC = gcc
CFLAGS = -std=c99 -c -Wall -Werror -fpic
LDFLAGS = -lm -lmesh -lfftw3 

all: numeric.so

numeric.so: libnumeric.o
	$(CC) -shared -o numeric.so libnumeric.o

libnumeric.o: libnumeric.c
	$(CC) $(CFLAGS) $(LDFLAGS) libnumeric.c

clean:
	rm -rf *o numeric.so

remake:
	make clean && make
