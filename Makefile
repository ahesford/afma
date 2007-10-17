CC=mpicc
FF=mpif90

OPTFLAGS?= -O3 -xP -parallel
CINCDIR= -I../../scaleme/include

CFLAGS?= $(OPTFLAGS) -Wall
FFLAGS?= $(OPTFLAGS) -Wall

OBJS= fsgreen.o integrate.o mlfma.o scaleme.o

all: $(OBJS)
	echo "Objects built."

clean:
	rm -f $(OBJS)

.SUFFIXES: .o .F .c

.F.o:
	$(FF) $(FFLAGS) -o $@ -c $<

.c.o:
	$(CC) $(CFLAGS) $(CINCDIR) -o $@ -c $<
