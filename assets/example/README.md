# Example data

A small **synthetic** DNA-encoded library in the [KinDEL](https://github.com/insitro/kindel) file layout,
for `-profile test` and the nf-test suites. It contains no real measurements. License: CC0-1.0.
Regenerate it with `uv run delt_pipeline simulate --outdir assets/example` (seeded, deterministic).

| File | Contents |
|---|---|
| `samplesheet.csv` | The three screens below with their metadata (`role`: `target` or `counter`) |
| `{mapk14,ddr1,bca}.parquet` | One selection each: 15 × 25 × 16 = 6,000 compounds, triplicate target (`seq_target_*`) and no-protein matrix (`seq_matrix_*`) counts, pre-selection counts (`seq_load`). Compounds with no counts are absent, as in KinDEL |
| `heldout/*_{on,off}dna.csv` | KD (nM) of 15 in-library compounds per kinase, plus 2 external compounds (no `molecule_hash`) |
| `synthons/reacted_{a,b,c}.csv` | Building blocks per position; the row number is the building-block id |

Planted signals (building-block ids, 0-based), see `src/delt_pipeline/library/simulate.py`:

- **MAPK14 binder family:** disynthon A2–B5, well synthesised, so it also has many control counts.
- **DDR1 binder family:** A4–B17. **BCA binder family:** A9–B20.
- **Kinase-shared family:** A7–B11 binds MAPK14 and DDR1 but not BCA.
- **Low-count outlier:** A1–B1–C1 has 6 MAPK14 target counts in one replicate and 0 control counts.
- **Bead binder:** A11–B3–C6 is enriched in the control and in every target.
- **Promiscuous compound:** A12–B8–C2 is enriched in all three proteins.

SMILES are simple fragments joined with `.`: valid, but not real DEL products.
