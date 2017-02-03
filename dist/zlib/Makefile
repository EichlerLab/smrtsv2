PWD=$(shell pwd)

all: lib/libz.so

lib/libz.so: zlib-1.2.8
	cd $< && ./configure --prefix=$(PWD) && make install clean

zlib-1.2.8: zlib-1.2.8.tar.gz
	tar zxvf $<

zlib-1.2.8.tar.gz:
	wget "https://github.com/madler/zlib/archive/v1.2.8.tar.gz" -O $@

.PHONY: clean

clean:
	rm -rf zlib-1.2.8*; rm -rf include lib share;