# Build SMRT-SV dependencies


# All
.PHONY: all
all: install_flags/dep_conda_link install_flags/dep_smrtsvtools_link bin/python


# python: Refuse unambiguous "python" interpreter (require "python2" or "python3")
bin/python: python_dummy/python_dummy.sh
	mkdir -p bin
	cp -f $< $@
	chmod 555 $@

# smrtsvtools
install_flags/dep_smrtsvtools_link: install_flags/dep_smrtsvtools_build
	mkdir -p install_flags
	ln -sf ../smrtsvtools/bin/alignfixup bin/alignfixup
	ln -sf ../smrtsvtools/bin/tlenadd bin/tlenadd
	date > $@

install_flags/dep_smrtsvtools_build: install_flags/dep_conda_build
	mkdir -p install_flags
	make -C smrtsvtools
	date > $@

# Conda
install_flags/dep_conda_link: install_flags/dep_conda_build
	mkdir -p bin install_flags

	# python3
	ln -sf ../conda/build/envs/python3/bin/python3 bin/python3
	ln -sf ../conda/build/envs/python3/bin/ipython3 bin/ipython3
	ln -sf ../conda/build/envs/python3/bin/snakemake bin/snakemake

	# python2
	ln -sf ../conda/build/envs/python2/bin/python2 bin/python2

	# pacbio
	ln -sf ../conda/build/envs/pacbio/bin/blasr bin/blasr
	ln -sf ../conda/build/envs/pacbio/bin/pbindex bin/pbindex
	ln -sf ../conda/build/envs/pacbio/bin/sawriter bin/sawriter
	ln -sf ../conda/build/envs/pacbio/bin/variantCaller bin/variantCaller
	ln -sf ../conda/build/envs/pacbio/bin/quiver bin/quiver
	ln -sf ../conda/build/envs/pacbio/bin/arrow bin/arrow

	# Tools
	ln -sf ../conda/build/envs/tools/bin/canu bin/canu
	ln -sf ../conda/build/envs/tools/bin/bedtools bin/bedtools
	ln -sf ../conda/build/envs/tools/bin/groupBy bin/groupBy
	ln -sf ../conda/build/envs/tools/bin/vcfkeepinfo bin/vcfkeepinfo
	ln -sf ../conda/build/envs/tools/bin/vcf2tsv bin/vcf2tsv
	ln -sf ../conda/build/envs/tools/bin/vcffixup bin/vcffixup
	ln -sf ../conda/build/envs/tools/bin/tabix bin/tabix
	ln -sf ../conda/build/envs/tools/bin/bgzip bin/bgzip
	ln -sf ../conda/build/envs/tools/bin/seqtk bin/seqtk
	ln -sf ../conda/build/envs/tools/bin/samtools bin/samtools
	ln -sf ../conda/build/envs/tools/bin/bamleftalign bin/bamleftalign
	#ln -sf ../conda/build/envs/tools/bin/RepeatMasker bin/RepeatMasker
	ln -sf ../conda/build/envs/tools/bin/blastp bin/blastp
	ln -sf ../conda/build/envs/tools/bin/trf bin/trf
	ln -sf ../conda/build/envs/tools/bin/bcftools bin/bcftools

	ln -sf ../conda/build/envs/tools/share/bwakit-0.7.15-1/bwa bin/bwa
	ln -sf ../conda/build/envs/tools/share/bwakit-0.7.15-1/k8 bin/k8
	ln -sf ../conda/build/envs/tools/share/bwakit-0.7.15-1/bwa-postalt.js bin/bwa-postalt.js

	#mkdir -p share
	#ln -sf ../conda/build/envs/tools/share/RepeatMasker share/

	# Flag complete
	date > $@

install_flags/dep_conda_build:
	mkdir -p install_flags
	make -C conda
	date > $@


# Clean
.PHONY: clean
clean:
	make -C conda clean
	make -C smrtsvtools clean
	rm -rf bin
	rm -rf install_flags
