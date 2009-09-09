LIBS += -framework Accelerate
OPTFLAGS= -O3 -march=core2 -mtune=core2 -D_MACOSX

universal: adbim afma

adbim afma:
	@echo "Building programs for $(ARCH)."
	make clean
	make -j4 $@32
	make clean
	make -j4 $@64
	@echo "Combining 32- and 64-bit $@."
	lipo -create $@64 $@32 -output $@
