PWD=$(shell pwd)

wgs-8.3rc2/Linux-amd64/bin: wgs-8.3rc2
	cd $</kmer && make LIBRARY_PATH=$(PWD)/../bz2/lib install && cd ../src && make LIBRARY_PATH=$(PWD)/../bz2/lib
	rm -f $@/blasr

wgs-8.3rc2: wgs-8.3rc2.tar.bz2
	tar jxvf $<

wgs-8.3rc2.tar.bz2:
	wget http://sourceforge.net/projects/wgs-assembler/files/wgs-assembler/wgs-8.3/wgs-8.3rc2.tar.bz2/download -O $@