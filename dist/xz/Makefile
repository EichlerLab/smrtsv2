PREFIX = $(abspath ..)

${PREFIX}/lib/liblzma.so: xz-5.2.3/configure
	cd $(dir $<) && ./configure --prefix=${PREFIX} && make install

xz-5.2.3/configure: xz-5.2.3.tar.gz
	tar -zxf $<
	touch $@

xz-5.2.3.tar.gz:
	wget https://tukaani.org/xz/xz-5.2.3.tar.gz
