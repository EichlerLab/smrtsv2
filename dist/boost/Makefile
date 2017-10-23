PREFIX=$(abspath ..)

${PREFIX}/lib/libboost_program_options.so: boost_1_65_1/b2
	cd boost_1_65_1 && ./b2 install

boost_1_65_1/b2: boost_1_65_1/bootstrap.sh
	cd boost_1_65_1 && ./bootstrap.sh --prefix=${PREFIX}

boost_1_65_1/bootstrap.sh: boost_1_65_1.tar.gz
	tar -zxf $<
	touch $@

boost_1_65_1.tar.gz:
	wget https://sourceforge.net/projects/boost/files/boost/1.65.1/boost_1_65_1.tar.gz
