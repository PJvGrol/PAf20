CC= bspcc
CFLAGS= -std=c99 -Wall -O3
LFLAGS= -lm

OBJIP= bspinprod.o bspedupack.o
OBJHK= bsphk.o bspedupack.o
OBJSSV= seqsieve.o bspedupack.o
OBJBSV= bspsieve.o bspedupack.o
OBJTW= bsptwins.o bspedupack.o
OBJCJ= bspconj.o bspedupack.o
OBJBEN= bspbench.o bspedupack.o
OBJSORT= bspsort_test.o bspsort.o bspedupack.o
OBJLU= bsplu_test.o bsplu.o bspedupack.o
OBJFFT= bspfft_test.o bspfft.o bspedupack.o
OBJMV= bspmv_test.o bspmv.o bspsparse_input.o bspedupack.o
OBJMATCH= bspmatch_test.o bspmatch.o bspsparse_input.o bspedupack.o

all: hk inprod ssieve sieve twins bench sort lu fft matvec match

hk: $(OBJHK)
	$(CC) $(CFLAGS) -o hk $(OBJHK) $(LFLAGS)

inprod: $(OBJIP)
	$(CC) $(CFLAGS) -o inprod $(OBJIP) $(LFLAGS)

ssieve: $(OBJSSV)
	$(CC) $(CFLAGS) -o ssieve $(OBJSSV) $(LFLAGS)

sieve: $(OBJBSV)
	$(CC) $(CFLAGS) -o sieve $(OBJBSV) $(LFLAGS)

twins: $(OBJTW)
	$(CC) $(CFLAGS) -o twins $(OBJTW) $(LFLAGS)

conj: $(OBJCJ)
	$(CC) $(CFLAGS) -o conj $(OBJCJ) $(LFLAGS)

bench: $(OBJBEN)
	$(CC) $(CFLAGS) -o bench $(OBJBEN) $(LFLAGS)

sort: $(OBJSORT)
	$(CC) $(CFLAGS) -o sort $(OBJSORT) $(LFLAGS)

lu: $(OBJLU)
	$(CC) $(CFLAGS) -o lu $(OBJLU) $(LFLAGS)

fft: $(OBJFFT)
	$(CC) $(CFLAGS) -o fft $(OBJFFT) $(LFLAGS)

matvec: $(OBJMV)
	$(CC) $(CFLAGS) -o matvec $(OBJMV) $(LFLAGS)

match: $(OBJMATCH)
	$(CC) $(CFLAGS) -o match $(OBJMATCH) $(LFLAGS)

bsphk.o: bspedupack.h
bspinprod.o:  bspedupack.h
seqsieve.o:  bspedupack.h
bspsieve.o:  bspedupack.h
bsptwins.o:  bspedupack.h
bspconj.o:  bspedupack.h
bspbench.o:   bspedupack.h
bspsort.o:    bspedupack.h
bspsort_test.o:    bspedupack.h
bsplu.o:    bspedupack.h
bsplu_test.o:    bspedupack.h
bspfft.o:    bspedupack.h
bspfft_test.o:    bspedupack.h
bspmv.o:    bspedupack.h
bspmv_test.o:    bspedupack.h bspsparse_input.h
bspmatch_test.o:    bspedupack.h bspsparse_input.h
bspsparse_input.o: bspedupack.h bspsparse_input.h
bspedupack.o: bspedupack.h

.PHONY: clean
clean:
	rm -f $(OBJHK) $(OBJIP) $(OBJSSV) $(OBJBSV) $(OBJTW) $(OBJCJ) $(OBJBEN) $(OBJSORT) $(OBJLU) $(OBJFFT) $(OBJMV) $(OBJMATCH) hk inprod ssieve sieve twins conj bench sort lu fft matvec match
