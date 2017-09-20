PREFIX = $(abspath ..)

${PREFIX}/lib/libhdf5_cpp.so: hdf5-1.8.19/configure
	cd $(dir $<) && ./configure --enable-cxx --prefix=${PREFIX} && make install

hdf5-1.8.19/configure: hdf5-1.8.19.tar.gz
	tar -zxvf $<
	touch $@

hdf5-1.8.19.tar.gz:
	wget https://support.hdfgroup.org/ftp/HDF5/current18/src/hdf5-1.8.19.tar.gz
