CC=mpicc
FF=mpif90
LD=$(FF)

LIBS= -lScaleME2 -lgmres -lfftw3f

OPTFLAGS= -fopenmp -O3 -march=core2 -mtune=core2
ARCHFLAGS= -arch x86_64 -arch i386 -D_MACOSX
ARCHLIBS= -framework Accelerate

CINCDIR= -I../scaleme/include -I/usr/local/include -I/opt/local/include
LIBDIR= -L../scaleme/fma2 -L../gmres -L/usr/local/lib -L/opt/local/lib

DFLAGS=
CFLAGS= $(OPTFLAGS) $(ARCHFLAGS) -Wall
FFLAGS= $(OPTFLAGS) $(ARCHFLAGS) -Wall
LFLAGS= $(OPTFLAGS) $(ARCHFLAGS)

OBJS= fsgreen.o integrate.o mlfma.o itsolver.o direct.o \
      excite.o io.o measure.o frechet.o cg.o

all: adbim afma tissue
	@echo "Building for Darwin (universal)."

afma: main.o $(OBJS)
	@echo "Building $@."
	$(LD) $(DFLAGS) $(LFLAGS) -o $@ $^ $(LIBDIR) $(LIBS) $(ARCHLIBS)

adbim: dbim.o $(OBJS)
	@echo "Building $@."
	$(LD) $(DFLAGS) $(LFLAGS) -o $@ $^ $(LIBDIR) $(LIBS) $(ARCHLIBS)

tissue: tissue.o
	@echo "Building $@."
	$(LD) $(DFLAGS) $(LFLAGS) -o $@ $^ $(LIBDIR) -lm

darwin32: OPTFLAGS= -fopenmp -O3 -march=prescott -mtune=prescott
darwin32: ARCHFLAGS= -arch i386 -D_MACOSX
darwin32: afma adbim tissue
	@echo "Building for Darwin (32 bits)."

bsd: OPTFLAGS= -fopenmp -O3 -mtune=native -march=native
bsd: ARCHLIBS= -lalapack_r -lptf77blas -lptcblas -latlas_r
bsd: ARCHFLAGS= -m64 -D_FREEBSD
bsd: afma adbim tissue
	@echo "Building for FreeBSD."

bluehive: OPTFLAGS= -O3 -openmp -xHost -parallel -no-prec-div
bluehive: ARCHLIBS= -lgslcblas -llapack -lblas -nofor_main
bluehive: ARCHFLAGS= -m64 -D_LINUX -I/usr/local/gsl/1.12-gnu4.1/include/gsl \
	-L/usr/local/gsl/1.12-gnu4.1/lib
bluehive: afma adbim tissue
	@echo "Building for Linux."

clean:
	rm -f $(OBJS) main.o dbim.o *.core core tissue.o

distclean: clean
	rm -f afma adbim tissue

.SUFFIXES: .o .f .c

.f.o:
	$(FF) $(DFLAGS) $(FFLAGS) -o $@ -c $<

.c.o:
	$(CC) $(DFLAGS) $(CFLAGS) $(CINCDIR) -o $@ -c $<
