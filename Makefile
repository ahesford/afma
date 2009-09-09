CC=mpicc
FF=mpif90
LD=mpif90

SCALEME_ROOT= $(HOME)/software/scaleme
LIBS= -lScaleME2 -lgmres

include mk/$(shell uname -s).mk

CINCDIR= -I$(SCALEME_ROOT)/include
LIBDIR= -L$(SCALEME_ROOT)/fma2 -L../gmres

CFLAGS= $(OPTFLAGS) -Wall
FFLAGS= $(OPTFLAGS) -Wall
LFLAGS= $(OPTFLAGS)

OBJS= fsgreen.o integrate.o mlfma.o scaleme.o itsolver.o \
	excite.o io.o measure.o frechet.o cg.o

all: adbim64 afma64

afma64 adbim64: override CFLAGS += -m64
afma64 adbim64: override FFLAGS += -m64
afma64 adbim64: override LFLAGS += -m64

afma32 adbim32: override CFLAGS += -m32
afma32 adbim32: override FFLAGS += -m32
afma32 adbim32: override LFLAGS += -m32

afma64 afma32: main.o $(OBJS)
	@echo "Building afma64."
	$(LD) $(LFLAGS) -o $@ $^ $(LIBDIR) $(LIBS)

adbim64 adbim32: dbim.o $(OBJS)
	@echo "Building adbim64."
	$(LD) $(LFLAGS) -o $@ $^ $(LIBDIR) $(LIBS)

clean:
	rm -f $(OBJS) main.o dbim.o *.core core

distclean: clean
	rm -f afma afma32 afma64 adbim adbim32 adbim64

.SUFFIXES: .o .f .c

.f.o:
	$(FF) $(FFLAGS) -o $@ -c $<

.c.o:
	$(CC) $(CFLAGS) $(CINCDIR) -o $@ -c $<
