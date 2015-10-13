export SMRT_ROOT=/net/eichler/vol24/projects/sequencing/pacbio/software/smrtanalysis
export PATH=`pwd`/bin:$PATH

# Use modules to load dependencies if that environment is available.
module list &> /dev/null

if [[ "$?" -eq "0" ]]
then
    module load anaconda/2.3.0
    module load samtools/1.2-complete

    module load hdf5/1.8.13
    module load netcdf/4.3.2
    module load R/3.1.0
    module load perl/5.14.2
    module load RepeatMasker/3.3.0

    module load mpfr/3.1.0 mpc/0.8.2 gmp/5.0.2 gcc/4.9.1
else
    export ANACONDA_ROOT=$HOME/anaconda
    export PATH=$ANACONDA_ROOT/envs/python2/bin:$ANACONDA_ROOT/envs/python3/bin:$PATH
fi
