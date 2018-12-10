# Setup conda environment: python3

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

pip install --upgrade pip

conda install -y \
    numpy==1.13.1 \
    pandas==0.20.3 \
    scipy==0.19.1 \
    pysam==0.15.1 \
    snakemake==5.3.0 \
    biopython==1.72 \
    ipython==7.2.0 \
    drmaa==0.7.9 \
    scikit-learn=0.19.0 \
    intervaltree==2.1.0
