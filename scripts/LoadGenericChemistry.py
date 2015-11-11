#!/usr/bin/env python

USAGE = \
"""
loadChemistry.py

 Load chemistry info into a cmp.h5 from a single bas.h5
 representative of all bas.h5 contributing to a local assembly.
 usage:
  % loadChemistry input.bax.h5 aligned_reads.cmp.h5
"""

import sys, h5py, numpy as np
from pbcore.io import *

class ChemistryLoadingException(BaseException): pass

STRING_DTYPE = h5py.special_dtype(vlen=bytes)

def safeDelete(group, dsName):
    if dsName in group:
        del group[dsName]

def writeTriples(movieInfoGroup, triple):
    movieNamesInCmpH5 = list(movieInfoGroup["Name"])

    safeDelete(movieInfoGroup, "BindingKit")
    safeDelete(movieInfoGroup, "SequencingKit")
    safeDelete(movieInfoGroup, "SoftwareVersion")

    shape = movieInfoGroup["Name"].shape
    bindingKit      = movieInfoGroup.create_dataset("BindingKit"     , shape=shape, dtype=STRING_DTYPE, maxshape=(None,))
    sequencingKit   = movieInfoGroup.create_dataset("SequencingKit"  , shape=shape, dtype=STRING_DTYPE, maxshape=(None,))
    softwareVersion = movieInfoGroup.create_dataset("SoftwareVersion", shape=shape, dtype=STRING_DTYPE, maxshape=(None,))

    for movieName in movieNamesInCmpH5:
        idx = movieNamesInCmpH5.index(movieName)
        bindingKit[idx]      = triple[0]
        sequencingKit[idx]   = triple[1]
        softwareVersion[idx] = triple[2]

    assert all(bindingKit.value      != "")
    assert all(sequencingKit.value   != "")
    assert all(softwareVersion.value != "")


def main():
    if len(sys.argv) != 3:
        print USAGE
        return -1

    baseFileName = sys.argv[1]
    cmpFname = sys.argv[2]

    f = h5py.File(cmpFname, "r+")
    movieInfoGroup = f["MovieInfo"]

    bas = BasH5Reader(baseFileName)
    chemTriple = bas.chemistryBarcodeTriple
    writeTriples(movieInfoGroup, chemTriple)

if __name__ == '__main__':
    main()
