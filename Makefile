CC=mpicc
FF=mpif90
LD=mpif90

SCALEME_ROOT= $(HOME)/software/scaleme

OPTFLAGS= -O3 -xP -parallel
CINCDIR= -I$(SCALEME_ROOT)/include
LIBDIR= -L$(SCALEME_ROOT)/fma2 -L/usr/local/lib

LIBS= -lScaleME2 -lgmres -llapack-icc -lblas-icc -lm

CFLAGS= $(OPTFLAGS) -Wall
FFLAGS= $(OPTFLAGS) -Wall
LFLAGS= $(OPTFLAGS) -nofor_main

OBJS= fsgreen.o integrate.o mlfma.o scaleme.o itsolver.o excite.o io.o \
      measure.o frechet.o cg.o

all: adbim afma
	@echo "Combined build."

adbim: $(OBJS) dbim.o
	@echo "Building acoustic DBIM."
	$(LD) $(LFLAGS) -o $@ dbim.o $(OBJS) $(LIBDIR) $(LIBS)

afma: $(OBJS) main.o
	@echo "Building acoustic MLFMA."
	$(LD) $(LFLAGS) -o $@ main.o $(OBJS) $(LIBDIR) $(LIBS)

clean:
	rm -f $(OBJS) main.o dbim.o adbim afma *.core core

.SUFFIXES: .o .f .c

.f.o:
	$(FF) $(FFLAGS) -o $@ -c $<

.c.o:
	$(CC) $(CFLAGS) $(CINCDIR) -o $@ -c $<
