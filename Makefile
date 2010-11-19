CC=mpicc
FF=mpif90
LD=$(FF)

LIBS= -lScaleME2 -lfftw3f

OPTFLAGS= -fopenmp -O3 -Xarch_x86_64 -march=core2 -Xarch_i386 -march=prescott
ARCHFLAGS= -arch x86_64 -arch i386 -D_MACOSX
ARCHLIBS= -framework Accelerate

CINCDIR= -I../scaleme/include -I/usr/local/include -I/opt/local/include
LIBDIR= -L../scaleme/fma2 -L/usr/local/lib -L/opt/local/lib

DFLAGS=
CFLAGS= $(OPTFLAGS) $(ARCHFLAGS)
FFLAGS= $(OPTFLAGS) $(ARCHFLAGS)
LFLAGS= $(OPTFLAGS) $(ARCHFLAGS)

FWDOBJS= main.o
INVOBJS= frechet.o cg.o dbim.o
OBJS= fsgreen.o integrate.o mlfma.o itsolver.o direct.o io.o measure.o util.o config.o

all: adbim afma tissue
	@echo "Building for Darwin (universal)."

afma: $(FWDOBJS) $(OBJS)
	@echo "Building $@."
	$(LD) $(DFLAGS) $(LFLAGS) -o $@ $^ $(LIBDIR) $(LIBS) $(ARCHLIBS)

adbim: $(INVOBJS) $(OBJS)
	@echo "Building $@."
	$(LD) $(DFLAGS) $(LFLAGS) -o $@ $^ $(LIBDIR) $(LIBS) $(ARCHLIBS)

tissue: tissue.o
	@echo "Building $@."
	$(LD) $(DFLAGS) $(LFLAGS) -o $@ $^ $(LIBDIR) -lm $(ARCHLIBS)

darwin32: OPTFLAGS= -fopenmp -O3 -march=prescott
darwin32: ARCHFLAGS= -arch i386 -D_MACOSX
darwin32: afma adbim tissue
	@echo "Building for Darwin (32 bits)."

bsd: OPTFLAGS= -fopenmp -O3 -mtune=native -march=native
bsd: ARCHLIBS= -lalapack_r -lptf77blas -lptcblas -latlas_r
bsd: ARCHFLAGS= -m64 -D_FREEBSD
bsd: afma adbim tissue
	@echo "Building for FreeBSD."

bluehive: CC= mpiicc
bluehive: FF= mpiifort
bluehive: OPTFLAGS= -O3 -openmp -xHost -parallel
bluehive: ARCHLIBS= -lgslcblas -llapack -lblas -nofor_main
bluehive: ARCHFLAGS= -m64 -D_LINUX -I/usr/local/gsl/1.12-gnu4.1/include/gsl \
	-L/usr/local/gsl/1.12-gnu4.1/lib
bluehive: afma adbim tissue
	@echo "Building for BlueHive Linux."

ranger: OPTFLAGS= -fastsse -mp
ranger: ARCHFLAGS= -D_LINUX -I$(TACC_GSL_INC)/gsl -I$(TACC_FFTW3_INC) \
	-L$(TACC_GSL_LIB) -L$(TACC_FFTW3_LIB)
ranger: ARCHLIBS= -Mnomain -lgslcblas -llapack -lblas
ranger: afma adbim tissue
	@echo "Building for TACC Ranger."

clean:
	rm -f $(OBJS) $(FWDOBJS) $(INVOBJS) *.core core tissue.o

distclean: clean
	rm -f afma adbim tissue

.SUFFIXES: .o .f .c

.f.o:
	$(FF) $(DFLAGS) $(FFLAGS) -o $@ -c $<

.c.o:
	$(CC) $(DFLAGS) $(CFLAGS) $(CINCDIR) -o $@ -c $<
