PREFIX = $(abspath ..)

${PREFIX}/bin/bwa: bwa-0.7.15/bwa
	mkdir -p ${PREFIX}/bin
	install $< ${PREFIX}/bin

bwa-0.7.15/bwa: bwa-0.7.15/Makefile
	cd bwa-0.7.15 && make

bwa-0.7.15/Makefile: bwa-0.7.15.tar.bz2
	tar -xjf bwa-0.7.15.tar.bz2

bwa-0.7.15.tar.bz2:
	wget https://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.15.tar.bz2
