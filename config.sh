export INSTALL_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${INSTALL_DIR}/dist/hdf5/lib

export PATH=${INSTALL_DIR}/bin:$PATH
export PATH=${INSTALL_DIR}/dist/miniconda/envs/python2/bin:${INSTALL_DIR}/dist/miniconda/envs/python3/bin:${INSTALL_DIR}/dist/miniconda/bin:${INSTALL_DIR}/dist/celera/wgs-8.3rc2/Linux-amd64/bin/:${INSTALL_DIR}/dist/amos-3.1.0/bin:${INSTALL_DIR}/canu/Linux-amd64/bin:$PATH
export PATH=/net/eichler/vol4/home/jlhudd/src/snakemake/bin:$PATH
