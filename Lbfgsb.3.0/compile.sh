#gcc -c routines.c; gcc routines.o -lm -lf2c -o DRIVER1 driver1.c


## changed timer code and added new function to calculate machine-precision 
gcc -c blas.c; gcc -c lbfgsb.c; gcc -c linpack.c; gcc -c timer.c
gcc timer.o blas.o linpack.o lbfgsb.o  -lm  -o DRIVER1 driver1.c


#gcc timer.o blas.o linpack.o lbfgsb.o  -lm -lf2c -o DRIVER1 driver1.c

