PREFIX=$(abspath ..)

${PREFIX}/bin/swig: swig-3.0.12/Makefile
	mkdir -p ${dir $@}
	cd swig-3.0.12 && make
	cd swig-3.0.12 && make install
	cd swig-3.0.12 && make clean

swig-3.0.12/Makefile: swig-3.0.12/configure
	cd swig-3.0.12 && ./configure --prefix=${PREFIX}

swig-3.0.12/configure: swig-3.0.12.tar.gz
	tar -zxvf swig-3.0.12.tar.gz
	touch $@

swig-3.0.12.tar.gz:
	wget http://downloads.sourceforge.net/project/swig/swig/swig-3.0.12/swig-3.0.12.tar.gz
