
.PHONY: install

install: Miniconda3-latest-Linux-x86_64.sh
	bash Miniconda3-latest-Linux-x86_64.sh -b -p $(shell pwd)
	./bin/conda create --yes --name python3 python=3
	source ./bin/activate python3
	pip install git+https://bitbucket.org/johanneskoester/snakemake.git
	pip install drmaa
	rm -f ./envs/python3/lib/python3.5/site-packages/snakemake/jobscript.sh
	ln -s $(shell pwd)/../snakemake/cluster/UGE/jobscript.sh ./envs/python3/lib/python3.5/site-packages/snakemake/jobscript.sh
	source bin/deactivate
	./bin/conda create --yes --name python2 python=2 biopython h5py networkx
	source ./bin/activate python2
	pip install intervaltree
	pip install pybedtools
	source bin/deactivate

Miniconda3-latest-Linux-x86_64.sh:
	wget https://repo.continuum.io/miniconda/$@