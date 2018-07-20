PREFIX = $(abspath ..)

${PREFIX}/lib/libz.so: zlib-1.2.11/configure
	cd $(dir $<) && ./configure --prefix=${PREFIX} && make install

zlib-1.2.11/configure: v1.2.11.tar.gz
	tar -zxf $<
	touch $@

v1.2.11.tar.gz:
	wget https://github.com/madler/zlib/archive/v1.2.11.tar.gz
