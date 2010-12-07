CC=mpicc
LD=$(CC)

LIBS= -lScaleME2 -lfftw3f -lm

OPTFLAGS= -fopenmp -O3 -Xarch_x86_64 -march=core2 -Xarch_i386 -march=prescott
ARCHFLAGS= -arch x86_64 -arch i386 -D_MACOSX
ARCHLIBS= -framework Accelerate

CINCDIR= -I../scaleme/include -I/usr/local/include -I/opt/local/include
LIBDIR= -L../scaleme -L/usr/local/lib -L/opt/local/lib

DFLAGS=
CFLAGS= $(OPTFLAGS) $(ARCHFLAGS)
LFLAGS= $(OPTFLAGS) $(ARCHFLAGS)

FWDOBJS= main.o
INVOBJS= frechet.o cg.o dbim.o
OBJS= fsgreen.o integrate.o mlfma.o itsolver.o direct.o io.o measure.o util.o config.o

EXECS= adbim afma tissue mat2grp

all: $(EXECS)
	@echo "Building for Darwin (universal)."

afma: $(FWDOBJS) $(OBJS)
	@echo "Building $@."
	$(LD) $(DFLAGS) $(LFLAGS) -o $@ $^ $(LIBDIR) $(LIBS) $(ARCHLIBS)

adbim: $(INVOBJS) $(OBJS)
	@echo "Building $@."
	$(LD) $(DFLAGS) $(LFLAGS) -o $@ $^ $(LIBDIR) $(LIBS) $(ARCHLIBS)

tissue: tissue.o
	@echo "Building $@."
	$(LD) $(DFLAGS) $(LFLAGS) -o $@ $^ $(LIBDIR) $(ARCHLIBS)

mat2grp: mat2grp.o
	@echo "Building $@."
	$(LD) $(DFLAGS) $(LFLAGS) -o $@ $^ $(LIBDIR) $(ARCHLIBS)

darwin32: OPTFLAGS= -fopenmp -O3 -march=prescott
darwin32: ARCHFLAGS= -arch i386 -D_MACOSX
darwin32: all
	@echo "Building for Darwin (32 bits)."

bsd: OPTFLAGS= -fopenmp -O3 -mtune=native -march=native
bsd: ARCHLIBS= -lalapack_r -lptf77blas -lptcblas -latlas_r
bsd: ARCHFLAGS= -m64 -D_FREEBSD
bsd: all
	@echo "Building for FreeBSD."

bluehive: CC= mpiicc
bluehive: OPTFLAGS= -O3 -openmp -xHost -parallel
bluehive: ARCHLIBS= -lgslcblas -llapack -lblas
bluehive: ARCHFLAGS= -m64 -D_LINUX -I/usr/local/gsl/1.12-gnu4.1/include/gsl \
	-L/usr/local/gsl/1.12-gnu4.1/lib
bluehive: all
	@echo "Building for BlueHive Linux."

ranger: OPTFLAGS= -fastsse -mp
ranger: ARCHFLAGS= -D_LINUX -I$(TACC_GSL_INC)/gsl -I$(TACC_FFTW3_INC) \
	-L$(TACC_GSL_LIB) -L$(TACC_FFTW3_LIB)
ranger: ARCHLIBS= -lgslcblas -llapack -lblas
ranger: all
	@echo "Building for TACC Ranger."

kraken: CC= cc
kraken: OPTFLAGS= -fastsse -mp
kraken: ARCHFLAGS= -D_LINUX -I$(GSL_DIR)/include/gsl -I$(FFTW_INC) \
	-L$(GSL_DIR)/lib -L$(FFTW_DIR)
kraken: ARCHLIBS= -lgslcblas -llapack -lblas
kraken: all
	@echo "Building for NICS Kraken."

clean:
	rm -f $(OBJS) $(FWDOBJS) $(INVOBJS) $(EXECS) tissue.o mat2grp.o
	rm -f *.core core 

.SUFFIXES: .o .c

.c.o:
	$(CC) $(DFLAGS) $(CFLAGS) $(CINCDIR) -o $@ -c $<
