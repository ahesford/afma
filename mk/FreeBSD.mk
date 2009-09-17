LIBS += -lalapack_r -lptf77blas -lptcblas -latlas_r
OPTFLAGS= -O3 -mtune=native -march=native -fopenmp
CINCDIR += -I/usr/local/lib
LIBDIR += -L/usr/local/lib
