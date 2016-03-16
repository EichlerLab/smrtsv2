export PATH=/usr/local/bin:/usr/bin:/bin
export INSTALL_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${INSTALL_DIR}/dist/hdf5/lib

# Use modules to load dependencies if that environment is available.
module list &> /dev/null

if [[ "$?" -eq "0" ]]
then
    module purge
    module load modules modules-init modules-gs/prod modules-eichler/prod
    #module load anaconda/2.3.0
    #module load samtools/1.2-complete
    module load bedops/2.4.0

    #module load hdf5/1.8.13
    #module load netcdf/4.3.2
    module load R/2.15.0
    module load perl/5.14.2
    module load RepeatMasker/3.3.0

    #module load mpfr/3.1.0 mpc/0.8.2 gmp/5.0.2 gcc/4.9.1
fi

export PATH=${INSTALL_DIR}/bin:$PATH
export PATH=${INSTALL_DIR}/dist/miniconda/envs/python2/bin:${INSTALL_DIR}/dist/miniconda/envs/python3/bin:${INSTALL_DIR}/dist/miniconda/bin:$PATH
