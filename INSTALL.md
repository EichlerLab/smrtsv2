# Install SMRT SV and dependencies

SMRT SV depends on Python 2.7 and Python 3 which are both available through the
free [Anaconda Scientific Python
Distribution](https://store.continuum.io/cshop/anaconda/).

## Get the code

Check out the repository with a recursive clone to fetch submodules.

```bash
mkdir ~/src/smrtsv
cd ~/src/smrtsv
git clone --recursive git@github.com:EichlerLab/pacbio_variant_caller.git .
```

## Install Python

Install the miniconda distribution of python.

```bash
cd dist/miniconda
sh install.sh
```

## Install all other dependencies

```bash
cd ~/src/smrtsv
make
```

## Install BLASR

### Install HDF5 with C++ support

Download latest [HDF5](http://www.hdfgroup.org) source code.

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

## Install PBcR

Download the [latest release of PBcR](http://wgs-assembler.sourceforge.net/wiki/index.php/PBcR).

```bash
wget http://sourceforge.net/projects/wgs-assembler/files/wgs-assembler/wgs-8.3/wgs-8.3rc2-Linux_amd64.tar.bz2/download -O wgs-8.3rc2-Linux_amd64.tar.bz2
```

Unpack binaries, create a central directory for Celera, and copy binaries
(except BLASR) into that directory.

```bash
tar jxvf wgs-8.3rc2-Linux_amd64.tar.bz2
cd wgs-8.3rc2/Linux-amd64/bin
rm -f blasr
sudo mkdir -p /usr/local/celera
rsync -arv * /usr/local/celera/
```

Update the `celera_dir` parameter in the `config.json` for your analysis to
point to the above directory.

## Install bedtools

Download the [latest release of bedtools](https://github.com/arq5x/bedtools2/releases/latest).

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

## Install samtools

Download the [latest release of samtools](https://github.com/samtools/samtools/releases/latest).

```bash
wget https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2
```

Unpack, build samtools, and install binaries, libraries, and headers into
central directories.

```bash
tar jxvf samtools-1.2.tar.bz2
cd samtools-1.2
make
sudo make install install-htslib
sudo rsync -arv *.h /usr/local/include/
```

## Install freebayes

Download the [freebayes source code](https://github.com/ekg/freebayes).

```bash
git clone --recursive git://github.com/ekg/freebayes.git
```

Select the latest version and build the binaries.

```bash
cd freebayes
git checkout v0.9.21 && git submodule update --recursive
make
```

## Install RepeatMasker

### Install sequence search engine (HMMER)

Download the [latest version of HMMER](http://hmmer.janelia.org/).

```bash
wget http://selab.janelia.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz
```

Unpack and build binaries.

```bash
tar zxvf hmmer-3.1b2-linux-intel-x86_64.tar.gz
cd hmmer-3.1b2-linux-intel-x86_64
./configure
make
sudo make install
```

### Install Tandem Repeat Finder (TRF)

Download the [latest version of TRF](http://tandem.bu.edu/trf/trf.download.html)
and copy the binary into a central directory on your PATH (e.g.,
/usr/local/bin).

### Install RepeatMasker

Download the [latest version of RepeatMasker](http://www.repeatmasker.org).

```bash
wget http://www.repeatmasker.org/RepeatMasker-open-4-0-5.tar.gz
```

Unpack the RepeatMasker scripts and libraries into a central directory.

```bash
sudo tar zxvf RepeatMasker-open-4-0-5.tar.gz --directory /usr/local
cd /usr/local/RepeatMasker
```

Configure RepeatMasker.

```bash
perl ./configure
```

Optionally, install [complete RepBase libraries](http://www.girinst.org/)
(registration required) as described in the [RepeatMasker installation
documentation](http://www.repeatmasker.org/RMDownload.html).

Finally, link the RepeatMasker script into a central binary path.

```bash
cd /usr/local/bin
sudo ln -s ../RepeatMasker/RepeatMasker .
```
