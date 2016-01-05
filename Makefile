############################
#  A makefile for dist  #
############################

SAMTOOLS  := $(shell samtools --version 2>/dev/null)
BEDTOOLS  := $(shell bedtools --version 2>/dev/null)
FREEBAYES := $(shell freebayes --version 2>/dev/null)
BLASR     := $(shell bin/blasr 2> /dev/null)
CELERA    := $(shell bin/PBcR 2> /dev/null)
JAVA      := $(shell bin/java 2> /dev/null)
QUIVER    := $(shell quiver --version 2> /dev/null)
PBH5TOOLS := $(shell cmph5tools.py --help 2> /dev/null)
PWD  = $(shell pwd)

all: checkBedtools checkSamtools checkFreebayes checkBlasr checkCelera checkJava checkQuiver checkPbH5Tools

#
# Install core genomics tools.
#

bedtools2:
	git submodule update --init dist/bedtools
	-cd dist/bedtools && make && make install prefix=$(PWD)

freebayes:
	git submodule update --init  dist/freebayes
	cd dist/freebayes && git submodule update --init
	cd dist/freebayes/vcflib && git submodule update --init
	-cd dist/freebayes && make
	-@ln -s ../dist/freebayes/bin/freebayes bin/freebayes
	-@ln -s ../dist/freebayes/bin/bamleftalign bin/bamleftalign

samtools:
	git submodule update --init dist/samtools
	git submodule update --init dist/htslib
	-cd dist/samtools && make
	-@ln -s ../dist/samtools/samtools bin/samtools

#
# Install Celera/PBcR and its dependencies.
#

dist/java:
	-cd $@ && make
	-@ln -s ../$@/jre1.8.0_65/bin/java bin/java

dist/celera:
	-cd $@ && make
	-@ln -s ../$@/wgs-8.3rc2/Linux-amd64/bin/PBcR bin/PBcR

#
# Install BLASR and its dependencies.
#

dist/hdf5:
	cd $@ && $(MAKE)

dist/zlib:
	cd $@ && $(MAKE)

dist/blasr: dist/hdf5 dist/zlib
	git submodule update --init $@
	-cd $@ && make HDF5INCLUDEDIR=$(PWD)/$</include HDF5LIBDIR=$(PWD)/$</lib LIBRARY_PATH=$(PWD)/$(word 2,$^)/lib:$(LIBRARY_PATH) && make install PREFIX=$(PWD) && make clean

#
# Install Quiver and its dependencies.
#

dist/swig:
	cd $@ && $(MAKE)

pbcore:
	git submodule update --init dist/$@
	-cd dist/$@ && source $(PWD)/dist/miniconda/bin/activate python2 && python setup.py install

ConsensusCore: dist/swig
	git submodule update --init dist/$@
	-cd dist/$@ && source $(PWD)/dist/miniconda/bin/activate python2 && python setup.py install --swig=$(PWD)/$</bin/swig --swig-lib=$(PWD)/$</share/swig/3.0.8

GenomicConsensus:
	git submodule update --init dist/$@
	-cd dist/$@ && source $(PWD)/dist/miniconda/bin/activate python2 && python setup.py install

#
# pbh5tools
#

pbh5tools:
	git submodule update --init dist/$@
	-cd dist/$@ && source $(PWD)/dist/miniconda/bin/activate python2 && python setup.py install

#
# Check for existing system-wide installations before building locally.
#

checkPbH5Tools:
ifdef PBH5TOOLS
	@echo "Found pbh5tools"
else
	@echo "Trying to install pbh5tools"
	$(MAKE) pbh5tools
endif

checkQuiver:
ifdef QUIVER
	@echo "Found Quiver version: $(QUIVER)"
else
	@echo "Trying to install Quiver"
	$(MAKE) pbcore ConsensusCore GenomicConsensus
endif

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

checkCelera:
ifdef CELERA
	@echo "Found Celera"
else
	@echo "Trying to install Celera"
	$(MAKE) dist/celera
endif

checkJava:
ifdef JAVA
	@echo "Found Java"
else
	@echo "Trying to install Java"
	$(MAKE) dist/java
endif
