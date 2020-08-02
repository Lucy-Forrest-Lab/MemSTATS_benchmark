# Benchmarking CE-Symm 2.0, SymD 1.61, AnAnaS 0.8, MSSD and QuatSymm
---

To benchmark the algorithms, download the MemSTATS dataset and pdb files, the scripts in this folder, as well as
the raw output of the symmetry-detection algorithms stored on [Zenodo](https://doi.org/10.5281/zenodo.3228539), and using Python 3.5, run:

    python benchmarking_with_MemSTATS.py

Results are written in the **results/** subfolder.

To use the tests stored in the **tests/** subfolder, run:

    pytest -v

from this folder. Note that the tests also require the pdb files, as well as the raw output of the symmetry-detection algorithms unpacked in this folder.
