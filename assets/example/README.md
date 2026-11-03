# Example data

Small **synthetic** datasets for `-profile test` and the nf-test suites. They contain no real
measurements. License: CC0-1.0.

| File | Contents |
|---|---|
| `samplesheet.csv` | The two datasets below, with their metadata |
| `compounds.csv.gz` | 40 compounds: SMILES, six descriptors (MolWt, LogP, TPSA, HBD, HBA, RotB) and a binary `active` label, with a few missing values. Gzipped |
| `compounds_dictionary.csv` | Data dictionary for `compounds`, with units (UO) and terms (CHEMINF) as IRIs |
| `expression.csv` | 30 samples × 12 genes (log2 expression), `batch` and `condition`. It has Windows line endings, a constant column (`spike_in`) and a mostly-missing column (`LINC00473`) to exercise cleaning |

Samplesheet paths are relative to the project root, so launch from there.
