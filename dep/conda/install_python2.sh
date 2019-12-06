# Setup conda environment: python2

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# Upgrading pip causes some issues (GitHub issue #3)
# pip install --upgrade pip

conda install -y \
    numpy=1.15.4 \
    scipy=1.1.0 \
    pandas=0.23.4 \
    pysam=0.15.1 \
    biopython=1.72 \
    intervaltree=2.1.0 \
    networkx=2.2 \
    pybedtools=0.8.0
