LIBS += -framework Accelerate
OPTFLAGS= -O3 -march=core2 -mtune=core2 -D_MACOSX

universal: adbim afma

adbim afma:
	make clean
	make -j4 $@32
	make clean
	make -j4 $@64
	lipo -create $@64 $@32 -output $@
