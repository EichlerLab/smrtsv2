#######################################
# Build SMRT-SV software dependencies #
#######################################

#
# Set environment
#

MAKEFILE_DIR = $(realpath $(dir $(firstword $(MAKEFILE_LIST))))

ifeq (${PATH},)
	PATH := ${MAKEFILE_DIR}/bin
else
	PATH := ${MAKEFILE_DIR}/bin:${PATH}
endif

ifeq (${LD_LIBRARY_PATH},)
	LD_LIBRARY_PATH := ${MAKEFILE_DIR}/dist/lib
else
	LD_LIBRARY_PATH := ${MAKEFILE_DIR}/dist/lib:${LD_LIBRARY_PATH}
endif

export PATH
export LD_LIBRARY_PATH


#
# Default rule: Bulid all
#

BIN_OUT = $(addprefix dist/bin/, samtools bcftools bwa bwa-postalt.js k8 seqtk samblaster bedtools vcffixup snakemake canu blasr variantCaller)

all: ${BIN_OUT}


#
# Library and bin dependencies
#

# htslib
dist/lib/libhts.so: dist/lib/liblzma.so
	cd dist/htslib && make

# zlib
dist/lib/libz.so:
	cd dist/zlib && make

# hdf5
dist/lib/libhdf5_cpp.so:
	cd dist/hdf5 && make

# xz (liblzma)
dist/lib/liblzma.so:
	cd dist/xz && make

# Boost
dist/lib/libboost_program_options.so: dist/bin/snakemake
	cd dist/boost && make

# Swig
dist/bin/swig:
	cd dist/swig && make


#
# SAM/BCF tools
#

# samtools
dist/bin/samtools: dist/lib/libhts.so
	cd dist/samtools && make

# bcftools
dist/bin/bcftools: dist/lib/libhts.so
	cd dist/bcftools && make


#
# Core genomics tools
#

# bwa
dist/bin/bwa:
	cd dist/bwa && make

# bwakit (seqtk, k8, and bwa-postalt.js)
dist/bin/bwa-postalt.js dist/bin/k8 dist/bin/seqtk dist/bin/samblaster:
	cd dist/bwakit && make

# bedtools
dist/bin/bedtools:
	cd dist/bedtools && make

# vcflib
dist/bin/vcffixup:
	cd dist/vcflib && make

# canu
dist/bin/canu:
	cd dist/canu && make


#
# Python environment with snakemake
#

# Python (2/3) and snakemake
dist/bin/snakemake:
	cd dist/miniconda && make


#
# PacBio toolkit
#

# BLASR (blasr sawriter samtoh5 samFilter)
dist/bin/blasr: dist/lib/libhdf5_cpp.so dist/lib/libz.so dist/bin/snakemake
	cd dist/blasr && make

# pbcore
dist/miniconda/envs/python2/lib/python2.7/site-packages/pbcore-1.4.0-py2.7.egg/pbcore/__init__.py: dist/bin/snakemake
	cd dist/pbcore && make

# ConsensusCore
dist/miniconda/envs/python2/lib/python2.7/site-packages/ConsensusCore-1.0.2-py2.7.egg/ConsensusCore.py: dist/miniconda/envs/python2/lib/python2.7/site-packages/pbcore-1.4.0-py2.7.egg/pbcore/__init__.py dist/lib/libboost_program_options.so dist/bin/swig
	cd dist/ConsensusCore && make

# GenomicConsensus (arrow and quiver) - For assembly polishing
dist/bin/variantCaller: dist/miniconda/envs/python2/lib/python2.7/site-packages/ConsensusCore-1.0.2-py2.7.egg/ConsensusCore.py
	cd dist/GenomicConsensus && make


# PRIV: Removed RepeatMasker for the private version (using module RepeatMasker/3.3.0)
#bin/RepeatMasker: bin/phmmer
#	-cd dist/RepeatMasker && make
#	-@ln -s ../dist/RepeatMasker/RepeatMasker/RepeatMasker bin/RepeatMasker
