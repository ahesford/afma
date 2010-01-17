CC=mpicc
FF=mpif90
LD=$(CC)

LIBS= -lScaleME2 -lgmres -lgfortran -lfftw3f -lm

OPTFLAGS= -fopenmp -O3 -march=core2 -mtune=core2
ARCHFLAGS= -arch x86_64 -arch i386 -D_MACOSX
ARCHLIBS= -framework Accelerate

CINCDIR= -I../scaleme/include -I/usr/local/include -I/opt/local/include
LIBDIR= -L../scaleme/fma2 -L../gmres -L/usr/local/lib -L/opt/local/lib

DFLAGS=
CFLAGS= $(OPTFLAGS) $(ARCHFLAGS) -Wall
FFLAGS= $(OPTFLAGS) $(ARCHFLAGS) -Wall
LFLAGS= $(OPTFLAGS) $(ARCHFLAGS)

OBJS= fsgreen.o integrate.o mlfma.o itsolver.o \
      excite.o io.o measure.o frechet.o cg.o

all: adbim afma
	@echo "Building for Darwin (universal)."

afma: main.o $(OBJS)
	@echo "Building $@."
	$(LD) $(DFLAGS) $(LFLAGS) -o $@ $^ $(LIBDIR) $(LIBS) $(ARCHLIBS)

adbim: dbim.o $(OBJS)
	@echo "Building $@."
	$(LD) $(DFLAGS) $(LFLAGS) -o $@ $^ $(LIBDIR) $(LIBS) $(ARCHLIBS)

darwin32: OPTFLAGS= -fopenmp -O3 -march=prescott -mtune=prescott
darwin32: ARCHFLAGS= -arch i386 -D_MACOSX
darwin32: afma adbim
	@echo "Building for Darwin (32 bits)."

bsd: OPTFLAGS= -fopenmp -O3 -mtune=native -march=native
bsd: ARCHLIBS= -lalapack_r -lptf77blas -lptcblas -latlas_r
bsd: ARCHFLAGS= -m64 -D_FREEBSD
bsd: afma adbim
	@echo "Building for FreeBSD."

linux: GSLDIR= /usr/local/gsl_gnu/lib
linux: ARCHLIBS= $(GSLDIR)/libgslcblas.a /usr/lib64/liblapack.so.3 /usr/lib64/libblas.so.3
linux: ARCHFLAGS= -m64 -D_LINUX
linux: CC= /usr/local/openmpi/bin/mpicc
linux: FF= /usr/local/openmpi/bin/mpif90
linux: afma adbim
	@echo "Building for Linux."

clean:
	rm -f $(OBJS) main.o dbim.o *.core core

distclean: clean
	rm -f afma adbim

.SUFFIXES: .o .f .c

.f.o:
	$(FF) $(DFLAGS) $(FFLAGS) -o $@ -c $<

.c.o:
	$(CC) $(DFLAGS) $(CFLAGS) $(CINCDIR) -o $@ -c $<
