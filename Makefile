CC=mpicc
FF=mpif90
LD=$(CC)

LIBS= -lScaleME2 -lgmres -lgfortran -lfftw3 -lm
CINCDIR= -I../scaleme/include -I/usr/local/include -I/opt/local/include
LIBDIR= -L../scaleme/fma2 -L../gmres -L/usr/local/lib -L/opt/local/lib

include mk/$(shell uname -s).mk

CFLAGS= $(OPTFLAGS) -Wall
FFLAGS= $(OPTFLAGS) -Wall
LFLAGS= $(OPTFLAGS)

OBJS= fsgreen.o integrate.o mlfma.o itsolver.o \
      excite.o io.o measure.o frechet.o cg.o

all: adbim afma

afma: main.o $(OBJS)
	@echo "Building $@."
	$(LD) $(LFLAGS) -o $@ $^ $(LIBDIR) $(LIBS)

adbim: dbim.o $(OBJS)
	@echo "Building $@."
	$(LD) $(LFLAGS) -o $@ $^ $(LIBDIR) $(LIBS)

clean:
	rm -f $(OBJS) main.o dbim.o *.core core

distclean: clean
	rm -f afma adbim

.SUFFIXES: .o .f .c

.f.o:
	$(FF) $(FFLAGS) -o $@ -c $<

.c.o:
	$(CC) $(CFLAGS) $(CINCDIR) -o $@ -c $<
