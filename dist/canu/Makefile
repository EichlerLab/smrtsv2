PREFIX = $(abspath ..)

${PREFIX}/bin/canu: canu-1.6/current/bin/canu
	ln -sf $(shell ../../scripts/relpath $(abspath $(dir $@)) $(abspath $<)) $@

canu-1.6/current/bin/canu: canu-1.6/src/Makefile
	mkdir -p current/bin
	ln -sf $(shell echo canu-1.6/*/bin/canu | head -n 1 | sed -re 's#canu-1.6/(.*)/bin/canu#\1#') canu-1.6/current
	make -C canu-1.6/src

canu-1.6/src/Makefile: v1.6.tar.gz
	tar -zxf v1.6.tar.gz
	touch $@

v1.6.tar.gz:
	wget https://github.com/marbl/canu/archive/v1.6.tar.gz
