# Setup conda environment: tools

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# Upgrading pip causes some issues (GitHub issue #3)
# pip install --upgrade pip

conda install -y \
    bcftools=1.9 \
    bedtools=2.29.0 \
    bwakit=0.7.15 \
    boost=1.70.0 \
    canu=1.8 \
    freebayes=1.3.1 \
    gxx_linux-64=7.3.0 \
    htslib=1.9 \
    seqtk=1.3 \
    samtools=1.9 \
    tabix=0.2.6 \
    vcflib=1.0.0_rc2