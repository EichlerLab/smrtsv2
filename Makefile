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
# Default rule: Build all
#

BIN_OUT = $(addprefix dist/bin/, samtools bcftools bwa bwa-postalt.js k8 seqtk samblaster bedtools vcffixup snakemake canu blasr variantCaller alignfixup minimap2)

all: ${BIN_OUT}


#
# Library and bin dependencies
#

# htslib
dist/lib/libhts.so: dist/lib/liblzma.so
	make -C dist/htslib

# zlib
dist/lib/libz.so:
	make -C dist/zlib

# hdf5
dist/lib/libhdf5_cpp.so:
	make -C dist/hdf5

# xz (liblzma)
dist/lib/liblzma.so:
	make -C dist/xz

# Boost
dist/lib/libboost_program_options.so: dist/bin/snakemake dist/lib/liblzma.so dist/lib/libz.so
	make -C dist/boost

# Swig
dist/bin/swig:
	make -C dist/swig


#
# SMRTSV Tools
#

# smrtsvtools
dist/bin/alignfixup: dist/lib/libboost_program_options.so dist/lib/libhts.so
	make -C dist/smrtsvtools

#
# SAM/BCF tools
#

# samtools
dist/bin/samtools: dist/lib/libhts.so
	make -C dist/samtools

# bcftools
dist/bin/bcftools: dist/lib/libhts.so
	make -C dist/bcftools


#
# Core genomics tools
#

dist/bin/minimap2:
	make -C dist/minimap2

# bwa
dist/bin/bwa:
	make -C dist/bwa

# bwakit (seqtk, k8, and bwa-postalt.js)
dist/bin/bwa-postalt.js dist/bin/k8 dist/bin/seqtk dist/bin/samblaster:
	make -C dist/bwakit

# bedtools
dist/bin/bedtools:
	make -C dist/bedtools

# vcflib
dist/bin/vcffixup:
	make -C dist/vcflib

# canu
dist/bin/canu:
	make -C dist/canu


#
# Python environment with snakemake
#

# Python (2/3) and snakemake
dist/bin/snakemake: dist/lib/libz.so
	make -C dist/miniconda


#
# PacBio toolkit
#

# BLASR (blasr sawriter samtoh5 samFilter)
dist/bin/blasr: dist/lib/libhdf5_cpp.so dist/lib/libz.so dist/bin/snakemake
	make -C dist/blasr

# pbcore
dist/miniconda/envs/python2/lib/python2.7/site-packages/pbcore-1.4.0-py2.7.egg/pbcore/__init__.py: dist/bin/snakemake
	make -C dist/pbcore

# ConsensusCore
dist/miniconda/envs/python2/lib/python2.7/site-packages/ConsensusCore-1.0.2-py2.7.egg/ConsensusCore.py: dist/miniconda/envs/python2/lib/python2.7/site-packages/pbcore-1.4.0-py2.7.egg/pbcore/__init__.py dist/lib/libboost_program_options.so dist/bin/swig
	make -C dist/ConsensusCore

# GenomicConsensus (arrow and quiver) - For assembly polishing
dist/bin/variantCaller: dist/miniconda/envs/python2/lib/python2.7/site-packages/ConsensusCore-1.0.2-py2.7.egg/ConsensusCore.py
	make -C dist/GenomicConsensus


# PRIV: Removed RepeatMasker for the private version (using module RepeatMasker/3.3.0)
#bin/RepeatMasker: bin/phmmer
#	-make -C dist/RepeatMasker
#	-@ln -s ../dist/RepeatMasker/RepeatMasker/RepeatMasker bin/RepeatMasker
