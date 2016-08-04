# Install SMRT-SV and dependencies

SMRT-SV requires git, Python (2.6.6 or later) and Perl (5.10.1 or later) for
installation.

## Get the code

Check out the repository with a recursive clone to fetch submodules and build SMRT-SV dependencies.

```bash
mkdir -p ~/src/smrtsv
cd ~/src/smrtsv
git clone --recursive git@github.com:EichlerLab/pacbio_variant_caller.git .
make
```

## Test installation

Print SMRT-SV help to confirm installation.

```bash
./bin/smrtsv.py --help
```

Copy the entire repository to your desired installation directory (e.g., `/usr/local/smrtsv`) and add that directory to your path.

```bash
export PATH=/usr/local/smrtsv:$PATH
```

Alternately, run `smrtsv.py` directly from the installation directory.

```bash
/usr/local/smrtsv/bin/smrtsv.py --help
```
