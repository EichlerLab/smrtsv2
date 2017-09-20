# Packages: python3

. bin/activate python3
pip install --upgrade pip
pip install numpy==1.13.1
pip install pandas==0.20.3
pip install scipy==0.19.1
pip install pysam==0.12
pip install snakemake==4.0.0
pip install biopython==1.70
pip install matplotlib==2.0.2 
pip install matplotlib-venn==0.11.5
pip install readline==6.2.4.1
pip install ipython==6.1.0
pip install drmaa==0.7.8
pip install psutil==5.2.2
pip install scikit-learn==0.19.0
pip install PyYAML==3.12
pip install bx-python==0.7.3
pip install openpyxl==2.5.0a3
pip install h5py==2.7.1
. bin/deactivate


# Packages python2

. bin/activate python2
pip install --upgrade pip
pip install numpy==1.13.1
pip install pandas==0.20.3
pip install scipy==0.19.1
pip install pysam==0.11
pip install biopython==1.70
pip install matplotlib==2.0.2
pip install matplotlib-venn==0.11.5
pip install readline==6.2.4.1
pip install ipython==5.3.0
pip install drmaa==0.7.8
pip install psutil==5.2.2
pip install scikit-learn==0.19.0
pip install PyYAML==3.12
pip install bx-python==0.7.3
pip install openpyxl==2.5.0a3
pip install h5py==2.7.1
. bin/deactivate

touch packages.flag
