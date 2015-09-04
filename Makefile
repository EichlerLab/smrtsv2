############################
#  A makefile for dist  #
############################

SAMTOOLS  := $(shell samtools --version 2>/dev/null)
BEDTOOLS  := $(shell bedtools --version 2>/dev/null)
FREEBAYES := $(shell freebayes --version 2>/dev/null)
PWD  = $(shell pwd)

all: checkBedtools checkSamtools checkFreebayes

bedtools2:
	git submodule update --init  dist/bedtools
	cd $(PWD)/dist/bedtools && make
	sleep 1
	-ln -s dist/bedtools/bin/bedtools bin/bedtools
freebayes:
	git submodule update --init  dist/freebayes
	cd $(PWD)/dist/freebayes && make
	sleep 1
	-ln -s dist/freebayes/bin/freebayes bin/freebayes
samtools:
	git submodule update --init dist/samtools
	cd $(PWD)/dist/samtools && make	
	sleep 1
	-ln -s dist/samtools/samtools bin/samtools
checkSamtools:
ifdef SAMTOOLS
	@echo "Found samtools version: $(SAMTOOLS)"
	-ln -s $(shell which samtools) bin/samtools
else
	@echo "Trying to install samtools"
	$(MAKE) samtools
endif

checkBedtools:
ifdef BEDTOOLS
	@echo "Found bedtools version: $(BEDTOOLS)"
	-ln -s $(shell which bedtools) bin/bedtools
else
	@echo "Trying to install bedtools"
	$(MAKE) bedtools
endif

checkFreebayes:
ifdef FREEBAYES
	@echo "Found freebayes version: $(FREEBAYES)"
	-ln -s $(shell which freebayes) bin/freebayes
else
	@echo "Trying to install freebayes"
	$(MAKE) freebayes
endif