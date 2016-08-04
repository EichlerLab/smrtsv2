############################
#  A makefile for dist  #
############################

PWD  = $(shell pwd)

all: bin/bedtools bin/samtools bin/bcftools bin/freebayes bin/blasr bin/PBcR bin/java dist/miniconda/envs/python2/bin/quiver dist/miniconda/envs/python2/bin/cmph5tools.py bin/RepeatMasker

#
# Install core genomics tools.
#

bin/RepeatMasker: bin/phmmer
	-cd dist/RepeatMasker && make
	-@ln -s ../dist/RepeatMasker/RepeatMasker/RepeatMasker bin/RepeatMasker

bin/phmmer:
	-cd dist/hmmer && make
	-@cp -f dist/hmmer/bin/* bin/

bin/bedtools:
	git submodule update --init dist/bedtools
	-cd dist/bedtools && make && make install prefix=$(PWD)

bin/freebayes:
	git submodule update --init  dist/freebayes
	cd dist/freebayes && git submodule update --init
	cd dist/freebayes/vcflib && git submodule update --init
	-cd dist/freebayes && make
	-@ln -s ../dist/freebayes/bin/freebayes bin/freebayes
	-@ln -s ../dist/freebayes/bin/bamleftalign bin/bamleftalign

bin/bcftools: dist/htslib/libhts.a
	git submodule update --init dist/bcftools
	git submodule update --init dist/htslib
	-cd dist/bcftools && $(MAKE)
	-@ln -s ../dist/bcftools/bcftools bin/bcftools

bin/samtools: dist/htslib/libhts.a
	git submodule update --init dist/samtools
	git submodule update --init dist/htslib
	-cd dist/samtools && $(MAKE)
	-@ln -s ../dist/samtools/samtools bin/samtools

dist/htslib/libhts.a:
	git submodule update --init dist/htslib
	-cd dist/htslib && $(MAKE)
	-@ln -s ../dist/htslib/bgzip bin/bgzip
	-@ln -s ../dist/htslib/tabix bin/tabix

#
# Install Celera/PBcR and its dependencies.
#

bin/java:
	-cd dist/java && $(MAKE)
	-@ln -s ../dist/java/jre1.8.0_65/bin/java bin/java

bin/PBcR: dist/bz2/bin/bzip2
	-cd dist/celera && $(MAKE)
	-@ln -s ../dist/celera/wgs-8.3rc2/Linux-amd64/bin/PBcR bin/PBcR

dist/bz2/bin/bzip2:
	-cd dist/bz2 && $(MAKE)

#
# Install BLASR and its dependencies.
#

dist/hdf5/lib/libhdf5_cpp.so:
	cd dist/hdf5 && $(MAKE)

dist/zlib/lib/libz.so:
	cd dist/zlib && $(MAKE)

bin/blasr: dist/hdf5/lib/libhdf5_cpp.so dist/zlib/lib/libz.so
	git submodule update --init dist/blasr
	-cd dist/blasr && make HDF5INCLUDEDIR=$(PWD)/dist/hdf5/include HDF5LIBDIR=$(PWD)/dist/hdf5/lib LIBRARY_PATH=$(PWD)/dist/zlib/lib:$(LIBRARY_PATH) && make install PREFIX=$(PWD) && make clean

#
# Install Quiver and its dependencies.
#

dist/swig/bin/swig:
	cd dist/swig && $(MAKE)

dist/miniconda/envs/python2/lib/python2.7/site-packages/pbcore-1.0.0-py2.7.egg: dist/miniconda/bin/activate
	git submodule update --init dist/pbcore
	-cd dist/pbcore && source $(PWD)/dist/miniconda/bin/activate python2 && sed -i 's/pysam == 0.8.1/pysam >= 0.8.1/' setup.py && python setup.py install && make clean

dist/miniconda/envs/python2/lib/python2.7/site-packages/ConsensusCore-1.0.0-py2.7.egg: dist/swig dist/swig/bin/swig dist/miniconda/bin/activate
	git submodule update --init dist/ConsensusCore
	-cd dist/ConsensusCore && source $(PWD)/dist/miniconda/bin/activate python2 && python setup.py install --swig=$(PWD)/$</bin/swig --swig-lib=$(PWD)/$</share/swig/3.0.8 && make clean

dist/miniconda/envs/python2/bin/quiver: dist/miniconda/envs/python2/lib/python2.7/site-packages/pbcore-1.0.0-py2.7.egg dist/miniconda/envs/python2/lib/python2.7/site-packages/ConsensusCore-1.0.0-py2.7.egg dist/miniconda/bin/activate
	git submodule update --init dist/GenomicConsensus
	-cd dist/GenomicConsensus && source $(PWD)/dist/miniconda/bin/activate python2 && sed -i 's/pysam == 0.8.1/pysam >= 0.8.1/' setup.py && python setup.py install && make clean
	touch $@

#
# pbh5tools
#

dist/miniconda/envs/python2/bin/cmph5tools.py: dist/miniconda/bin/activate
	git submodule update --init dist/pbh5tools
	-cd dist/pbh5tools && source $(PWD)/dist/miniconda/bin/activate python2 && python setup.py install

#
# Install miniconda
#

dist/miniconda/bin/activate:
	cd dist/miniconda && sh install.sh
