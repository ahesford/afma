CC=mpicc
LD=$(CC)

FFTW= fftw3f
LIBS= -lScaleME -l$(FFTW)

OPTFLAGS= -fopenmp -O3 -march=core2
ARCHFLAGS= -D_MACOSX -flax-vector-conversions
ARCHLIBS= -framework Accelerate

CINCDIR= -I../scaleme/include -I/usr/local/include
LIBDIR= -L../scaleme -L/usr/local/lib

DFLAGS=
CFLAGS= $(OPTFLAGS) $(ARCHFLAGS)
LFLAGS= $(OPTFLAGS) $(ARCHFLAGS)

FWDOBJS= main.o
INVOBJS= frechet.o cg.o dbim.o
OBJS= fsgreen.o integrate.o mlfma.o itsolver.o direct.o io.o measure.o util.o config.o

EXECS= adbim afma tissue mat2grp lapden

all: $(EXECS)
	@echo "Building for Darwin (64-bit)."

afma: $(FWDOBJS) $(OBJS)
	@echo "Building $@."
	$(LD) $(DFLAGS) $(LFLAGS) -o $@ $^ $(LIBDIR) $(LIBS) $(ARCHLIBS)

adbim: $(INVOBJS) $(OBJS)
	@echo "Building $@."
	$(LD) $(DFLAGS) $(LFLAGS) -o $@ $^ $(LIBDIR) $(LIBS) $(ARCHLIBS)

tissue: tissue.o
	@echo "Building $@."
	$(LD) $(DFLAGS) $(LFLAGS) -o $@ $^ $(LIBDIR) $(ARCHLIBS)

lapden: lapden.o
	@echo "Building $@."
	$(LD) $(DFLAGS) $(LFLAGS) -o $@ $^ $(LIBDIR) \
		-l$(FFTW)_threads $(LIBS) $(ARCHLIBS)

mat2grp: mat2grp.o
	@echo "Building $@."
	$(LD) $(DFLAGS) $(LFLAGS) -o $@ $^ $(LIBDIR) $(ARCHLIBS)

bsd: LD= mpif77
bsd: OPTFLAGS= -fopenmp -O3 -mtune=native -march=native
bsd: ARCHLIBS= -lalapack_r -lptf77blas -lptcblas -latlas_r
bsd: ARCHFLAGS= -m64 -D_ATLAS
bsd: all
	@echo "Building for FreeBSD."

bluehive: CC= mpiicc
bluehive: LD= mpiifort
bluehive: GSL_DIR= /usr/local/gsl/1.12-gnu4.1
bluehive: OPTFLAGS= -O3 -openmp -xHost -parallel
bluehive: ARCHLIBS= -lgslcblas -llapack -lblas -nofor_main
bluehive: ARCHFLAGS= -m64 -D_LINUX -I$(GSL_DIR)/include -L$(GSL_DIR)/lib
bluehive: all
	@echo "Building for BlueHive Linux."

ranger: LD= mpif77
ranger: OPTFLAGS= -fastsse -mp
ranger: ARCHFLAGS= -D_LINUX -I$(TACC_GSL_INC)/gsl -I$(TACC_FFTW3_INC)
ranger: ARCHLIBS= -L$(TACC_GSL_LIB) -L$(TACC_FFTW3_LIB) \
	-Mnomain -lgslcblas -llapack -lblas
ranger: all
	@echo "Building for TACC Ranger."

kraken: CC= cc
kraken: LD= ftn
kraken: OPTFLAGS= -fastsse -mp
kraken: ARCHFLAGS= -D_LINUX -I$(GSL_DIR)/include/gsl -I$(FFTW_INC)
kraken: ARCHLIBS= -L$(GSL_DIR)/lib -L$(FFTW_DIR) -Mnomain -lgslcblas
kraken: all
	@echo "Building for NICS Kraken."

habis: LD= mpif77
habis: OPTFLAGS= -fopenmp -O3 -mtune=native -march=native
habis: ARCHLIBS= -L$(ATLAS_DIR)/lib -L$(FFTW_DIR)/lib \
	-llapack -lptf77blas -lptcblas -latlas
habis: ARCHFLAGS= -m64 -D_ATLAS -I$(ATLAS_DIR)/include -I$(FFTW_DIR)/include
habis: all
	@echo "Building for HABIS."

clean:
	rm -f $(OBJS) $(FWDOBJS) $(INVOBJS) $(EXECS) tissue.o mat2grp.o lapden.o
	rm -f *.core core 

.SUFFIXES: .o .c

.c.o:
	$(CC) $(DFLAGS) $(CFLAGS) $(CINCDIR) -o $@ -c $<
