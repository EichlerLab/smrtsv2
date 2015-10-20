PWD=$(shell pwd)

all: lib/libhdf5_cpp.so

lib/libhdf5_cpp.so:
	wget http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.15-patch1.tar.gz -O hdf5-1.8.15-patch1.tar.gz && \
	tar zxvf hdf5-1.8.15-patch1.tar.gz && \
	cd hdf5-1.8.15-patch1 && \
	./configure --enable-cxx --prefix=$(PWD) && $(MAKE) install clean \
	rm -rf hdf5-1.8.15-patch1*
