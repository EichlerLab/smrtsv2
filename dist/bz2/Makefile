PWD=$(shell pwd)

bin/bzip2 : bzip2-1.0.6
	cd $< && make && make install PREFIX=$(PWD) && make clean

bzip2-1.0.6 : bzip2-1.0.6.tar.gz
	tar zxvf $<

bzip2-1.0.6.tar.gz :
	wget http://www.bzip.org/1.0.6/$@ -O $@