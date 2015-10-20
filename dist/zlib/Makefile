PWD=$(shell pwd)

all: lib/libz.so

lib/libz.so:
	wget http://zlib.net/zlib-1.2.8.tar.gz && \
	tar zxvf zlib-1.2.8.tar.gz && \
	cd zlib-1.2.8 && \
	./configure --prefix=$(PWD) && \
	$(MAKE) install clean
	rm -rf zlib-1.2.8*
