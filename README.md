# MemSTATS_benchmark
---

MemSTATS (Membrane protein Structures And Their Symmetries) is a manually-curated benchmark set of **quaternary and internal symmetries** in integral membrane protein structures. A protein is included in the dataset only if its complex or one of its membrane-spanning chains has a **distinct structural fold** compared to all other proteins within the benchmark. Moreover, a symmetry description is included only if the structural repeats are at least partially **within the membrane-embedded region of the protein**. 

MemSTATS consists of an unambiguous and easily-parsable PDB file for each protein in the benchmark, a general table with all symmetry descriptions, a table with only internal (or intrasubunit) symmetry descriptions of unique chains, and a table with only quaternary symmetry descriptions of the complexes. Each entry provides the following general structural information:

•	“Fold-Abbreviation” – an abbreviation of the structural fold, of which the given protein is a representative

•	“Fold-Name” – name of the structural fold, of which the given protein is a representative

•	“PDB” – the PDB code of the protein structure

•	“Resolution” – the resolution of the structure

•	“Secondary-Structure-Composition” – whether the protein is primarily composed of transmembrane alpha helices or beta strands

•	“TM-Chains” – transmembrane chains within the structure

•	“Structurally-Unique-Chains” – the names of the protein transmembrane chains that have distinct folds compared to all the rest of the dataset

•	“Structure-Title” – name of the given protein

•	“Reference” – a number pointing to a reference that describes the structure or its symmetry (see "References.txt")


While compiling the set, we aimed to describe the diversity of membrane protein symmetries, while minimizing trivial cases. The following features are described for each symmetry (note that “;” is used to separate each distinct symmetry description within a structure): 

•	“Order” – the number of repeats that constitute a given symmetry

•	“Repeats-Topology” – the topology of the repeats with respect to the membrane 

•	“Closed/Open” – whether the symmetry can be described by a point group (closed, such as C2 or D4) or has a translational component (open, such as helical or linear symmetries)

•	“Internal/Quaternary” – whether the symmetry is intrasubunit (internal) or intersubunit (quaternary)

•	 “Described” – whether the symmetry is explicitly mentioned in the literature

•	“Repeats” – approximate amino-acid range of the repeats (should be considered accurate within 20 amino acids on the overall). Each bracket defines the regions that have been aligned to each other. To illustrate how the repeat ranges should be interpreted, consider the following example:

> (B_2-290,C_2-217)(E_2-290,F_2-217)(H_2-290,I_2-217)

This entry corresponds to two repeats associated with a single symmetry. Repeat 1 constitutes of *residues 2-290 in chain B* AND *residues 2-290 in chain E* AND *residues 2-290 in chain H*. Repeat 2 constitutes of *residues 2-217 in chain C* AND *residues 2-217 in chain F* AND *residues 2-217 in chain I*.

•	“Interdigitating” – whether the repeats are interdigitating 

The PDB files of the structures are included in MemSTATS_pdbs.tar.gz. The files are formatted so that:

•	the biological unit and the names of the chains are consistent with the benchmark dataset

•	there is only one 'altloc' allowed to avoid ambiguity in individual algorithms processing.

