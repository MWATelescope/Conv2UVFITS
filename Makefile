# Makefile for corr2uvfits.

CFLAGS=-g -O -Wall -D_FILE_OFFSET_BITS=64 -L. 
CFITSIO_INCS=$(shell pkg-config --silence-errors --cflags cfitsio)
CFITSIO_LIBS=$(shell pkg-config --silence-errors --libs cfitsio)

TARGETS=corr2uvfits Lphasor test_readuvfits

all: $(TARGETS)

corr2uvfits: corr2uvfits.c convutils.c convutils.h uvfits.c uvfits.h
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS) -o corr2uvfits corr2uvfits.c convutils.c uvfits.c $(CFITSIO_LIBS) -lcfitsio -lstarlink_pal -lm

Lphasor: Lphasor.c convutils.c convutils.h uvfits.c
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS) -o Lphasor Lphasor.c convutils.c uvfits.c -lstarlink_pal $(CFITSIO_LIBS) -lcfitsio -lm

test_readuvfits: test_readuvfits.c uvfits.c
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS) -o $@ $^ $(CFITSIO_LIBS) -lcfitsio -lstarlink_pal -lm

clean:
	rm -f *.o $(TARGETS) 
	rm -rf corr2uvfits.dSYM
