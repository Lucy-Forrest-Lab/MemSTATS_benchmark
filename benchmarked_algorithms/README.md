# Symmetry-detection algorithms tested on the structures in the MemSTATS dataset
---

The folder **symmetry/** contains the raw outputs of processing the pdb files in the MemSTATS dataset symmetry detection algorithms (SymD 1.61, CE-Symm 2.0, AnAnaS 0.8). In the case of QuatSymm (also referred to as BioJava/RCSB), the data has been manually extracted from the PDB database.

The Python dictionary stored in *.tm_archive.pkl* contains information about transmembrane regions of the MemSTATS structures obtained by processing them with [PPM](https://dx.doi.org/10.1093%2Fnar%2Fgkr703).

To obtain the structure files and the MemSTATS dataset, as well as the code for benchmarking the algorithms against MemSTATS, visit the [MemSTATS_benchmark GitHub](https://github.com/AntoniyaAleksandrova/MemSTATS_benchmark).

Note that to use the files without issues, your file system needs to be case-sensitive, i.e. it should be able to distinguish between a file name *FILE.pdb* and a file name *file.pdb*.
