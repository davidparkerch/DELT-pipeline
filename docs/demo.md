# Five-minute walkthrough

A script for showing this project live, and the notes behind each step. Every number comes from `just run-kindel`
on KinDEL and can be reproduced; the commands are in the README. Open the [explorer](https://davidparkerch.github.io/delt-pipeline/explorer.html)
and the [report](https://davidparkerch.github.io/delt-pipeline/) in two tabs first.

## 1. The problem (30 seconds)

A DEL selection returns sparse, noisy counts: most compounds have none, and the rest mix binding with synthesis
yield, PCR and sticking to the beads. The job is to turn them into a short list worth resynthesising. KinDEL is a
public 81 M-compound library with three selections (MAPK14, DDR1, carbonic anhydrase) and measured KD for about
40 resynthesised compounds per kinase. That makes it possible to check a method, not only to run it.

## 2. Does it work? The validation table (60 seconds)

*Show the table in the README.*

- The count baseline matches the paper exactly in all four columns; the Poisson score is within 0.02 in three of four.
  The fourth (MAPK14 off-DNA, 0.24 against 0.17) is a difference I report, not hide.
- With 30 to 49 compounds a correlation has a standard error near 0.15. So the five scores are **not** separable,
  and the pipeline reports all of them rather than naming a winner.
- The useful signal is elsewhere: scores that use the control (EB, ER, Poisson) beat the target-only z-score off-DNA,
  because the z-score fills MAPK14's top 2,000 with 1,744 bead or counter-screen binders.

## 3. Why one count is not evidence (45 seconds)

*Show the saturation figure.*

- MAPK14 is sampled at 0.82 counts per compound (65.8 M counts, the figure in the paper, over 80.7 M compounds). Protocols
  ask for ten.
- Half of the compounds seen at all were seen once. So the pipeline ranks by the **lower bound** of an exact interval,
  and a planted test shows a 6-count outlier beating a real family under naive enrichment but not under the bound.
- This is the RNA-seq library-complexity curve, and the same idea as a size factor, applied to DEL.

## 4. Hits are families (45 seconds)

*Show the cube in the explorer, then the replicate plot.*

- Replicates agree at r ≈ 0.00 to 0.08 per compound and r ≈ 0.86 per disynthon. Lines in the cube are disynthons.
- Triage asks for replicate support and an enriched disynthon scored **without the compound itself**, plus flags for
  the BCA counter-screen and the beads. Shared kinase hits are reported, not penalised: plausible hinge binders.

## 5. The explorer: mark, search, open (60 seconds)

1. Choose DDR1. Tick only tier 1.
2. Lasso the top right of the scatter plot. The ring appears in the cube and the check marks in the table: three views, one set of marked hits.
3. Search `S(=O)(=O)N` for **building blocks containing**. Nine blocks match. Point at the BCA column: the two para-sulfonamides have
   AUC 0.996 and 0.990 there and 0.44 to 0.59 in the kinases. Carbonic anhydrase is inhibited by primary sulfonamides, and the pipeline
   found them without being told (ranks 1 and 2 of 1,127; about 2 in a million by chance).
4. Click a row: structure, per-replicate counts, five scores, flags, reasons, and each building block with its set-enrichment AUC.
5. Export the marked hits to CSV.

## 6. From 2,000 hits to a plate (30 seconds)

Tick **Shortlist only**. The 96 picks come from 96 structural clusters. The top 96 of the ranking come from 17
(MAPK14) because an enriched family carries all its analogues up together. The cost is stated: median rank 441 against 50.

## 7. What a model adds, and what it does not (30 seconds)

On-DNA, the best statistic (0.78) beats the best fingerprint model (0.57 for MAPK14, 0.74 for DDR1). Off-DNA the model helps only
if it is given the control (0.56 against −0.13), and no more than the noise. Leave out every disynthon of the validated compounds and
MAPK14 on-DNA drops from 0.52 to 0.02. The paper saw the same pattern.

## 8. The DESeq2 analogy, and where it breaks (30 seconds)

Transfers: size factors become exposures; replicates estimate noise; `lfcShrink` becomes the empirical-Bayes shrunken ratio.
Breaks: median-of-ratios needs most features to be non-zero (here almost all counts are 0); three near-zero replicates cannot
estimate a per-compound dispersion; there are no gene-level priors, but the building-block structure is the real way to
share strength. The control is a bead experiment, not an untreated baseline.

## Questions to expect

- *Why Poisson and not negative binomial?* Overdispersion is real (see the replicate agreement). The pipeline handles it with a
  prior fitted on random building blocks and by asking for more than a bound of 1; a hierarchical model over building blocks
  would be the next step and is listed in the limitations.
- *Why not use the pre-selection counts as the normaliser?* KinDEL warns that its noise propagates for low-abundance compounds.
  The pipeline uses them as a yield proxy per building block, and the matrix control as the background.
- *What about truncated products and linker effects?* Not observable in these counts. The `low_load` flag is a proxy.
- *Is the KinDEL sub-library biased?* The pipeline keeps every building block of the validated compounds plus 60 random ones;
  library-wide priors are fitted on the random ones only.
