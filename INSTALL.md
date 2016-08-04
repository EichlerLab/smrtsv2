# Install SMRT SV and dependencies

SMRT SV depends on Python 2.7 and Python 3 which are both available through the
free [Anaconda Scientific Python
Distribution](https://store.continuum.io/cshop/anaconda/).

## Get the code

Check out the repository with a recursive clone to fetch submodules.

```bash
mkdir -p ~/src/smrtsv
cd ~/src/smrtsv
git clone --recursive git@github.com:EichlerLab/pacbio_variant_caller.git .
```

## Install with make

```bash
cd ~/src/smrtsv
make
```

## Install RepeatMasker

SMRT-SV expects RepeatMasker to be installed in your $PATH to annotate the
repeat content of SVs. Use the following instructions to install RepeatMasker
with HMMER for the alignment engine.

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
