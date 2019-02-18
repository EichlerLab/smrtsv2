# Setup conda environment: pacbio

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

pip install --upgrade pip

conda install -y \
    genomicconsensus==2.3.2 \
    blasr==5.3.2 \
    hdf5==1.10.2
