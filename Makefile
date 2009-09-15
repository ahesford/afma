CC=mpicc
FF=mpif90
LD=mpif90

LIBS= -lScaleME2 -lgmres
CINCDIR= -I../scaleme/include
LIBDIR= -L../scaleme/fma2 -L../gmres

include mk/$(shell uname -s).mk

CFLAGS= $(OPTFLAGS) -Wall
FFLAGS= $(OPTFLAGS) -Wall
LFLAGS= $(OPTFLAGS)

OBJS= fsgreen.o integrate.o mlfma.o parfmm.o itsolver.o \
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
