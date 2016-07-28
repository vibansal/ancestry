CC=gcc -D_GNU_SOURCE

all:
	$(MAKE) -C parsebam/samtools-0.1.18 all
	$(MAKE) -C parsebam all
	$(MAKE) -C . ancestry

ancestry: lbfgsb.o ancestry-gll.c Lbfgsb.3.0/lbfgsb.h Lbfgsb.3.0/lbfgsb.c
	$(CC) lbfgsb.o -lm -o ANCESTRY ancestry-gll.c 


lbfgsb.o: Lbfgsb.3.0/lbfgsb.h Lbfgsb.3.0/lbfgsb.c
	$(CC) -o lbfgsb.o -c Lbfgsb.3.0/lbfgsb.c

clean:
	rm -f lbfgsb.o ANCESTRY
	$(MAKE) -C parsebam clean
	#$(MAKE) -C parsebam/samtools-0.1.18 clean
