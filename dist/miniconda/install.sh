#! /bin/bash

wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -f -b -p `pwd`
./bin/conda create --yes --name python3 python=3 pandas pysam
#./bin/conda install -n python3 pip
. ./bin/activate python3
pip install --upgrade pip
pip install git+https://bitbucket.org/snakemake/snakemake.git
pip install drmaa
rm -f ./envs/python3/lib/python3.5/site-packages/snakemake/jobscript.sh
ln -s `pwd`/../snakemake/cluster/UGE/jobscript.sh ./envs/python3/lib/python3.5/site-packages/snakemake/jobscript.sh
. bin/deactivate
./bin/conda create --yes --name python2 python=2 biopython h5py networkx scipy
. ./bin/activate python2
pip install intervaltree
pip install pandas
pip install pysam==0.8.4
pip install pybedtools
. bin/deactivate
