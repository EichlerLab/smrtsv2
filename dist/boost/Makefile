PREFIX=$(abspath ..)

${PREFIX}/lib/libboost_program_options.so: boost_1_61_0/b2
	cd boost_1_61_0 && ./b2 install

boost_1_61_0/b2: boost_1_61_0/bootstrap.sh
	cd boost_1_61_0 && ./bootstrap.sh --prefix=${PREFIX}

boost_1_61_0/bootstrap.sh: boost_1_61_0.tar.gz
	tar -zxvf $<
	touch $@

boost_1_61_0.tar.gz:
	wget https://sourceforge.net/projects/boost/files/boost/1.61.0/boost_1_61_0.tar.gz
