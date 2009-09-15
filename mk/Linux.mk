GSLDIR= /usr/local/gsl_gnu/lib
LIBS += $(GSLDIR)/libgslcblas.a /usr/lib64/liblapack.so.3 /usr/lib64/libblas.so.3
OPTFLAGS= -O3 -march=nocona -mtune=nocona -m64
override CC= /usr/local/openmpi/bin/mpicc
override FF= /usr/local/openmpi/bin/mpif90
