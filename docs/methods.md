# Methods

How DELT-pipeline turns DEL selection counts into ranked, validated hits, and why each step is done
the way it is. Code references point to `src/delt_pipeline/`.

## 1. Data: KinDEL

[KinDEL](https://github.com/insitro/kindel) (Chen et al., ICML 2025) is a public trisynthon DEL:
377 × 1,127 × 190 building blocks (about 81 M compounds). It was screened against MAPK14 and DDR1, two kinases,
and BCA (bovine carbonic anhydrase), a non-kinase. Each selection was run in triplicate, alongside a
no-protein **matrix** control (beads only) in triplicate and a **pre-selection** ("load") sample of the
naive library. Each screen table has 68 M rows; compounds with no counts at all are absent.

**Counts are sequence counts, not read counts.** Every KinDEL compound is encoded by several redundant
DNA tags, and the counts are the number of *distinct* tags seen. PCR duplicates of one tag count once. This is
the DEL equivalent of UMI deduplication: MAPK14's 514 M reads collapse to 66 M sequence counts, about
8 PCR copies per molecule. The decoding step (split reads by selection index, match codons to building
blocks with error-tolerant barcode matching, collapse duplicates) happened upstream. This pipeline
starts from the decoded counts.

For each kinase, KinDEL also provides measured KD values for 30–50 compounds that were resynthesised **on-DNA**
(fluorescence polarisation) and **off-DNA** (SPR). These are the ground truth for [validation](#5-validation).

## 2. Sub-library (`library/subset.py`)

A full screen (68 M rows, about 7 GB of Parquet) is never loaded into memory. polars streams it, with a peak
of about 4 GB, so three screens fit on a 16 GB machine:

1. **Whole-library depth first.** The first pass totals every count column over all 68 M compounds.
   These totals are the exposures that normalise every rate. If they came from the subset instead,
   every rate would shift relative to the published analysis, and the validation would stop being
   like-for-like.
2. **Choose BB-As** (`select_sublibrary`). The subset keeps every BB-A used by a held-out compound,
   so all validation compounds are present, plus random BB-As up to `--n_bb_a` (default 60, about 11 M
   compounds per screen). Slicing by building block keeps the complete B × C structure, which synthon
   analysis needs; a random sample of compounds would destroy it.
3. **Map synthon SMILES to integer building-block ids** and write a compact, sorted Parquet table.
   KinDEL also lists compounds with an *empty* position (`BLANK` or an empty string): truncated
   products with only two building blocks. That is why the paper counts 378 × 1,128 × 191 while the synthon
   tables have 377 × 1,127 × 190. These compounds count towards the whole-library depth but are
   excluded from scoring and reported as `n_unmapped` (about 0.6 % of the MAPK14 sub-library).

**Bias, stated plainly:** BB-As chosen because a validated compound uses them are enriched for binders.
Per-compound scores do not depend on the rest of the library, but a prior fitted to "the library" does.
So library-wide quantities (the empirical-Bayes prior) are fitted **only on compounds whose BB-A was
drawn at random** (`random_bb_a` column).

## 3. Normalisation

- **Depth:** each arm's exposure is the sum of its three replicates' whole-library totals, E_t and E_c.
- **Why not DESeq2 size factors?** Median-of-ratios needs features that are non-zero in every sample.
  Here over 80 % of compounds have zero counts in any one replicate, so the median ratio is undefined
  or 0. Total counts are the only stable size factor.
- **Pre-selection load** is reported per compound, but it is itself a sparse Poisson sample of the naive library.
  Dividing by it adds noise to exactly the low-count compounds that are already noisy. Target and control are
  compared for the same compound, so load cancels out of the target-vs-control ratio.

## 4. Enrichment scores (`analyse/enrichment.py`)

With pooled counts c_t (target) and c_c (control) per compound:

| Score | Formula | Strength | Weakness |
|---|---|---|---|
| `count_diff` | mean(c_t) − mean(c_c) | KinDEL's "Counts" baseline; favours abundant binders | Ignores depth and noise |
| `er` | ((c_t + α)/E_t) / ((c_c + α)/E_c) | Intuitive, depth-normalised | **Ignores evidence**: 3 vs 0 beats 200 vs 40 |
| `poisson_lower` | Lower bound of the exact CI of the rate ratio | Needs both a high ratio and enough counts | Conservative for rare compounds |
| `z_faver` | (p_o − p_e)/√(p_e(1 − p_e)), p_e = 1/D | Depth-independent (Faver 2019) | Target-only: bead binders score high too |
| `eb_log2_ratio` | Posterior mean log2 ratio under a Beta prior | Shrinks low-count noise, keeps strong signals | Needs a sensible prior population |

**Exact rate-ratio interval.** Conditional on the total n = c_t + c_c, the target count is binomial:
c_t | n ~ Bin(n, π), with π = R·E_t / (R·E_t + E_c) and R the rate ratio. Clopper–Pearson bounds on π map
exactly to bounds on R:

    R = π / (1 − π) · E_c / E_t

"No enrichment" is π₀ = E_t / (E_t + E_c), **not** 0.5. With twice the target depth, 200 vs 100 counts
is *not* enrichment, and a test pins this (`test_poisson_interval_uses_the_exposure_not_one_half`).
Ranking by the lower bound is the Poisson-enrichment idea of Gerry et al. (2019), the "Poisson"
baseline in KinDEL.

**Empirical-Bayes shrinkage.** π gets a Beta(a, b) prior, fitted by maximising the Beta-Binomial likelihood
over the random-BB-A compounds. Identical (c_t, c_c) pairs are grouped, so millions of compounds reduce to
a few thousand likelihood terms. The posterior is Beta(a + c_t, b + c_c), and its mean log-odds has a
closed form, ψ(a + c_t) − ψ(b + c_c), so there is no sampling and no per-compound optimisation.

**The low-count outlier** (`tests/python/test_score.py::test_planted_signal`). The simulated library
plants a compound with 6 target counts in one replicate and none in the control. Naive ER ranks it
above the best-supported member of a real binder family; the Poisson lower bound and EB shrinkage rank
it below. This is the single most important behaviour of a DEL scorer.

## 5. Validation (`analyse/validate.py`)

For every score, we compute the **negative Spearman correlation** with measured KD (lower KD = tighter
binder, so higher is better). This is the metric of the KinDEL paper, so our "Counts" and "Poisson" columns
can be compared with its published baselines (MAPK14 on-DNA 0.778 / 0.737, DDR1 on-DNA 0.695 / 0.779).

Every measurement row counts, including the two DDR1 off-DNA compounds that were measured twice, exactly as
in the KinDEL benchmark code. With that, the "Counts" baseline reproduces all four published values exactly.

Held-out compounds that are absent from a screen table had no counts in any sample. They get the score
of an all-zero compound instead of being dropped; dropping them would remove the true negatives and
inflate every correlation.

## 6. DESeq2 thinking: what transfers and what breaks

DEL hit calling is count-based differential abundance: target vs control instead of treated vs untreated.

**What transfers**

- **Size factors → exposures.** Different sequencing depths must be divided out before comparing arms.
- **Replicates estimate noise.** Triplicates show whether a signal is reproducible (`target_support`,
  replicate correlations in the report). On MAPK14, target replicates correlate at r ≈ 0.01–0.08 per
  compound, but at r ≈ 0.85–0.87 once counts are summed per A–B disynthon. Replicates agree on
  synthons, not on single compounds.
- **Shrinkage of low-count fold changes.** `lfcShrink` exists because a fold change estimated from 3 reads
  is meaningless. `eb_log2_ratio` is the same idea with a Beta-Binomial prior.
- **MA-plot intuition.** Ratio against mean count shows that extreme ratios live at low counts.

**What breaks**

- **Extreme sparsity.** RNA-seq genes have hundreds of reads; most DEL compounds have 0. Median-of-ratios
  normalisation fails, and the noise model must handle mostly zeros.
- **Per-feature dispersion is not estimable.** DESeq2 shares information across genes to estimate each gene's
  dispersion. Three near-zero replicates per compound carry almost no information, so the prior has to
  be library-wide.
- **No gene-level priors, but a combinatorial structure.** Genes have annotations; compounds have building
  blocks. The real "sharing of strength" in DEL comes from synthons: a disynthon seen across 190 different
  third building blocks is far more convincing than one bright compound. That is why the cube plot and
  synthon aggregation matter more than any single-compound statistic.
- **The control is not a baseline condition.** The matrix control measures binding to the beads, not
  an untreated state. Compounds enriched in the control, or in an unrelated protein like BCA, are
  artefacts, not "down-regulated".

## 7. Limitations

- **On-DNA vs off-DNA.** DEL counts measure binding while attached to DNA. Off-DNA KD correlates only loosely,
  as the KinDEL paper also reports.
- **No competition arm.** KinDEL has no selection with a known inhibitor, so allosteric and orthosteric binders
  cannot be separated.
- **Sub-library.** Results cover about 16 % of BB-As. The pipeline runs on the full library with
  `--n_bb_a 377` given enough memory and time.
- **GPU.** Not used. Every step is vectorised NumPy/SciPy on CPU, and outputs must be byte-stable for snapshot
  tests. A structure-based model (fingerprints → Poisson loss) is a natural next step and could use a GPU.

## References

- Chen B., Danel T., et al. *KinDEL: DNA-Encoded Library Dataset for Kinase Inhibitors.* ICML 2025, arXiv:2410.08938.
- Gerry C. J., Wawer M. J., Clemons P. A., Schreiber S. L. *DNA barcoding a complete matrix of stereoisomeric small molecules.* JACS 141, 10225 (2019).
- Faver J. C., et al. *Quantitative comparison of enrichment from DNA-encoded chemical library selections.* ACS Comb. Sci. 21, 75 (2019).
- Love M. I., Huber W., Anders S. *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.* Genome Biol. 15, 550 (2014).
- Favalli N., Bassi G., Scheuermann J., Neri D. *DNA-encoded chemical libraries – achievements and remaining challenges.* FEBS Lett. 592, 2168 (2018).
