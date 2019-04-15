# Setup conda environment: pacbio

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# Upgrading pip causes some issues (GitHub issue #3)
# pip install --upgrade pip

conda install -y \
    bam2fastx==1.3.0 \
    bax2bam==0.0.9 \
    blasr==5.3.2 \
    genomicconsensus==2.3.2 \
    hdf5==1.10.2
