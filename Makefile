CC=mpicc
FF=mpif90
LD=mpif90

SCALEME_ROOT= ../../scaleme

OPTFLAGS?= -O3 -xP -parallel
CINCDIR= -I$(SCALEME_ROOT)/include
LIBDIR= -L$(SCALEME_ROOT)/fma2

LIBS?= -lScaleME2 -llapack -lblas -lm

CFLAGS?= $(OPTFLAGS) -Wall
FFLAGS?= $(OPTFLAGS) -Wall
LFLAGS?= $(OPTFLAGS)

OBJS= fsgreen.o integrate.o mlfma.o scaleme.o itsolver.o excite.o io.o \
      measure.o main.o cgmres.o

afma: $(OBJS)
	$(LD) $(LFLAGS) -o $@ $(OBJS) $(LIBDIR) $(LIBS)

clean:
	rm -f $(OBJS) afma *.core core

.SUFFIXES: .o .f .c

.f.o:
	$(FF) $(FFLAGS) -o $@ -c $<

.c.o:
	$(CC) $(CFLAGS) $(CINCDIR) -o $@ -c $<
