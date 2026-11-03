"""Column names of KinDEL tables and of the pipeline's compact count tables."""

POSITIONS = ("a", "b", "c")
"""Synthon positions (cycles) of the trisynthon library."""

REPLICATES = (1, 2, 3)

KINDEL_TARGET = tuple(f"seq_target_{rep}" for rep in REPLICATES)
KINDEL_CONTROL = tuple(f"seq_matrix_{rep}" for rep in REPLICATES)
KINDEL_LOAD = "seq_load"
KINDEL_COUNTS = (*KINDEL_TARGET, *KINDEL_CONTROL, KINDEL_LOAD)
"""Sequence-count columns of a KinDEL screen (distinct encodings per molecule, i.e. PCR-deduplicated)."""

TARGET = tuple(f"target_{rep}" for rep in REPLICATES)
CONTROL = tuple(f"matrix_{rep}" for rep in REPLICATES)
LOAD = "load"
COUNTS = (*TARGET, *CONTROL, LOAD)
RENAME = dict(zip(KINDEL_COUNTS, COUNTS, strict=True))
"""KinDEL count column -> compact name (``seq_target_1`` -> ``target_1``, ``seq_load`` -> ``load``)."""

BB_COLUMNS = tuple(f"bb_{position}" for position in POSITIONS)
SMILES_COLUMNS = tuple(f"smiles_{position}" for position in POSITIONS)
ID = "molecule_hash"
