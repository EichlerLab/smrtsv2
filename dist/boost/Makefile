PWD=$(shell pwd)

lib/libboost_system.a: boost_1_61_0
	cd $< && ./bootstrap.sh && ./b2 --prefix=$(PWD) install

boost_1_61_0: boost_1_61_0.tar.gz
	tar zxvf $<

boost_1_61_0.tar.gz:
	wget http://downloads.sourceforge.net/project/boost/boost/1.61.0/$@
