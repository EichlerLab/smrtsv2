
PWD=$(shell pwd)

all: bin/phmmer

bin/phmmer: hmmer-3.1b2
	cd $< && ./configure --prefix=$(PWD) && make && make install && make clean

hmmer-3.1b2: hmmer-3.1b2.tar.gz
	tar zxvf $<

hmmer-3.1b2.tar.gz:
	wget http://eddylab.org/software/hmmer3/3.1b2/$@ -O $@
