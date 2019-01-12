# Setup conda environment: tools

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

pip install --upgrade pip

conda install -y \
    canu==1.8 \
    bedtools==2.27.1 \
    boost==1.67.0 \
    htslib==1.9 \
    gcc=4.8.5 \
    vcflib==1.0.0_rc1 \
    tabix==0.2.6 \
    bwakit=0.7.15 \
    seqtk=1.3 \
    freebayes=1.2.0 \
    repeatmasker=4.0.7 \
    trf=4.09 \
    bcftools=1.9
