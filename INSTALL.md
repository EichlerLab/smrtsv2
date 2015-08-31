# Install SMRT SV and dependencies

SMRT SV depends on Python 2.7 and Python 3 which are both available through the
free [Anaconda Scientific Python
Distribution](https://store.continuum.io/cshop/anaconda/).

## Install Python

Download the latest Anaconda installation (e.g.,
"Anaconda-2.3.0-Linux-x86_64.sh" for 64-bit Linux systems). Run the installation
script and follow on-screen instructions.

```bash
bash Anaconda-2.3.0-Linux-x86_64.sh
```

Create a conda environment for Python 2.

```bash
conda create --name python2 biopython
```

Create a conda environment for Python 3.

```bash
conda create --name python3 python=3
```

Install [Snakemake](https://bitbucket.org/johanneskoester/snakemake/) in the
Python 3 environment.

```bash
source activate python3
pip install snakemake
source deactivate
```

Install [intervaltree](https://pypi.python.org/pypi/intervaltree/) in the Python
2 environment.

```bash
source activate python2
pip install intervaltree
source deactivate
```

Finally, prepare your environment to run Snakemake with Python 3 and all other
scripts with Python 2. This works by placing Python 2 binaries on the system
PATH before Python 3 binaries. The following example assumes Anaconda is
installed in the current user's home directory.

```bash
export PATH=$HOME/anaconda/envs/python2/bin:$HOME/anaconda/envs/python3/bin:$PATH
```

## Install bedtools

Download [latest release of bedtools](https://github.com/arq5x/bedtools2/releases/latest).

```bash
wget https://github.com/arq5x/bedtools2/releases/download/v2.24.0/bedtools-2.24.0.tar.gz
```

Unpack and build the binaries.

```bash
tar zxvf bedtools-2.24.0.tar.gz
cd bedtool2
make
sudo make install
```

## Install BLASR

### Install HDF5 with C++ support

Download latest HDF5 source code.

```bash
wget http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.15-patch1.tar.gz
```

Unpack source and build binaries and libraries with C++ support enabled.

```bash
tar zxvf hdf5-1.8.15-patch1.tar.gz
cd hdf5-1.8.15-patch1
./configure --enable-cxx --prefix=/usr/local/hdf5
make
sudo make install
```

### Install zlib

Download the latest source code for [zlib](http://www.zlib.net/).

```bash
wget http://zlib.net/zlib-1.2.8.tar.gz
```

Unpack and build zlib libraries.

```bash
tar zxvf zlib-1.2.8.tar.gz
cd zlib-1.2.8
./configure
make
sudo make install
```

### Build BLASR

Setup HDF5 include and library directories.

```bash
export HDF5INCLUDEDIR=/usr/local/hdf5/include
export HDF5LIBDIR=/usr/local/hdf5/lib
```

Download BLASR source code and build binaries.

```bash
git clone https://github.com/EichlerLab/blasr.git
cd blasr
make
sudo make install PREFIX=/usr/local
```
