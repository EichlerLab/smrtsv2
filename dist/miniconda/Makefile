PREFIX = $(abspath ..)

ALL_BIN = $(addprefix ${PREFIX}/bin/, python python3 ipython ipython3 snakemake python2 ipython2)


.PHONY: install
install: ${ALL_BIN} ${PREFIX}/bin_py2/python

${PREFIX}/bin_py2/python: packages.flag
	mkdir -p $(dir $@)
	ln -sf $(abspath .)/envs/python2/bin/python2 $@

${PREFIX}/bin/python: packages.flag
	mkdir -p $(dir $@)
	ln -sf $(abspath .)/envs/python3/bin/python3 ${PREFIX}/bin/python

${PREFIX}/bin/python3: packages.flag
	mkdir -p $(dir $@)
	ln -sf $(abspath .)/envs/python3/bin/python3 ${PREFIX}/bin/python3

${PREFIX}/bin/ipython3: packages.flag
	mkdir -p $(dir $@)
	ln -sf $(abspath .)/envs/python3/bin/ipython3 ${PREFIX}/bin/ipython3

${PREFIX}/bin/ipython: packages.flag
	mkdir -p $(dir $@)
	ln -sf $(abspath .)/envs/python3/bin/ipython3 ${PREFIX}/bin/ipython

${PREFIX}/bin/snakemake: packages.flag
	mkdir -p $(dir $@)
	ln -sf $(abspath .)/envs/python3/bin/snakemake ${PREFIX}/bin/snakemake

${PREFIX}/bin/python2: packages.flag
	mkdir -p $(dir $@)
	ln -sf $(abspath .)/envs/python2/bin/python2 ${PREFIX}/bin/python2

${PREFIX}/bin/ipython2: packages.flag
	mkdir -p $(dir $@)
	ln -sf $(abspath .)/envs/python2/bin/ipython2 ${PREFIX}/bin/ipython2


packages.flag: envs/python2/bin/python2.7 envs/python3/bin/python3.6
	sh install_packages.sh

envs/python2/bin/python2.7: bin/conda
	bin/conda create --yes --name python2 python=2.7.13
	touch $@

envs/python3/bin/python3.6: bin/conda
	bin/conda create --yes --name python3 python=3.6.2
	touch $@

bin/conda: Miniconda3-4.3.21-Linux-x86_64.sh
	bash Miniconda3-4.3.21-Linux-x86_64.sh -f -b -p $(realpath .)

Miniconda3-4.3.21-Linux-x86_64.sh:
	wget https://repo.continuum.io/miniconda/Miniconda3-4.3.21-Linux-x86_64.sh
