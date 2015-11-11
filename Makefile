############################
#  A makefile for dist  #
############################

SAMTOOLS  := $(shell samtools --version 2>/dev/null)
BEDTOOLS  := $(shell bedtools --version 2>/dev/null)
FREEBAYES := $(shell freebayes --version 2>/dev/null)
BLASR     := $(shell bin/blasr 2> /dev/null)
PWD  = $(shell pwd)

all: checkBedtools checkSamtools checkFreebayes checkBlasr

bedtools2:
	git submodule update --init  dist/bedtools
	-cd dist/bedtools && make && make install prefix=$(PWD)

freebayes:
	git submodule update --init  dist/freebayes
	cd dist/freebayes && git submodule update --init
	cd dist/freebayes/vcflib && git submodule update --init
	-cd dist/freebayes && make
	-@ln -s ../dist/freebayes/bin/freebayes bin/freebayes

samtools:
	git submodule update --init dist/samtools
	git submodule update --init dist/htslib
	-cd dist/samtools && make
	-@ln -s ../dist/samtools/samtools bin/samtools

dist/blasr: dist/hdf5 dist/zlib
	git submodule update --init $@
	-cd $@ && make HDF5INCLUDEDIR=$(PWD)/$</include HDF5LIBDIR=$(PWD)/$</lib LIBRARY_PATH=$(PWD)/$(word 2,$^)/lib:$(LIBRARY_PATH) && make install PREFIX=$(PWD) && make clean

dist/hdf5:
	cd $@ && $(MAKE)

dist/zlib:
	cd $@ && $(MAKE)

checkSamtools:
ifdef SAMTOOLS
	@echo "Found samtools version: $(SAMTOOLS)"
	-@ln -s $(shell which samtools) bin/samtools
else
	@echo "Trying to install samtools"
	$(MAKE) samtools
endif

checkBedtools:
ifdef BEDTOOLS
	@echo "Found bedtools version: $(BEDTOOLS)"
	-@ln -s $(shell which bedtools) bin/bedtools
else
	@echo "Trying to install bedtools"
	$(MAKE) bedtools2
endif

checkFreebayes:
ifdef FREEBAYES
	@echo "Found freebayes version: $(FREEBAYES)"
	-@ln -s $(shell which freebayes) bin/freebayes
else
	@echo "Trying to install freebayes"
	$(MAKE) freebayes
endif

checkBlasr:
ifdef BLASR
	@echo "Found BLASR"
else
	@echo "Trying to install BLASR"
	$(MAKE) dist/blasr
endif
