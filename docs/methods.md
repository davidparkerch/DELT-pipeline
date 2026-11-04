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

## 3a. Sampling depth and saturation (`library/depth.py`)

A selection can only say something about a compound that was sampled. The protocol literature asks for
sequence counts that **exceed the library size, ideally tenfold** (Decurtins et al., *Nat. Protoc.* 2016).
The same scan that totals the whole-library depth also builds, for every count column and for the pooled
target and matrix counts, a histogram of compounds per exact count. From it the pipeline reports (in
`<id>.counts_qc.json`, key `library.saturation`) the counts per library member, the fraction of the library
observed, the share of observed compounds seen exactly once, and a **rarefaction curve**.

The rarefaction curve is the DEL version of an RNA-seq library-complexity curve. Keeping each count with
probability *f* leaves a compound with *k* counts unobserved with probability (1 − *f*)^*k*, so the expected
number of observed compounds is Σ *n*_k (1 − (1 − *f*)^*k*). That is exact under binomial thinning and needs no
resampling. It cannot extrapolate beyond the real depth, so a dotted line shows what a perfectly uniform
library would do at higher depth.

| MAPK14 | Counts | Counts per library member | Library observed | Uniform library would show |
|---|---|---|---|---|
| Target, pooled | 65.8 M | **0.82** | 35.9 % | 55.7 % |
| Matrix control, pooled | 72.4 M | 0.90 | 35.7 % | 59.2 % |
| Pre-selection load | 563 M | 7.0 | 81.9 % | 99.9 % |

What this means:

- **The total of 65,767,904 target counts equals the 65.8 M sequence counts the KinDEL paper reports**
  (its Appendix A.2), spread over 80.7 M possible compounds. The target selection is sampled at 0.82 counts
  per compound, about twelve times below the recommended depth.
- **Half of the observed compounds were seen exactly once** (pooled target: 49.9 %). A single count carries almost
  no information about enrichment, which is why the pipeline ranks by a *lower confidence bound* and asks for
  disynthon support instead of ranking by counts. It also predicts the near-zero per-compound replicate
  correlation in section 4: each replicate is sampled at 0.1 to 0.4 counts per compound.
- **The library is not uniform.** Even the deeply sequenced pre-selection load, at 7 counts per member, sees only
  81.9 % of the library where a uniform library would show 99.9 %. Roughly one compound in five is missing or
  rare before any selection, because of synthesis yield, truncation or loss. A compound without counts in
  the target is therefore **not evidence of non-binding**.
- The matrix and load columns are identical in the three released screens (the same control experiment),
  so their rows repeat across screens.

The reads behind these sequence counts are not in the released tables, so the PCR duplication rate (reads per
distinct encoding, 6 to 8.5 in the paper) is not computed here.

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

## 6. Synthon aggregation (`analyse/synthons.py`)

Counts are summed per **monosynthon** (one building block, any partners) and per **disynthon** (a pair
of building blocks, any third), in all three orientations: AB, AC and BC. Each pooled group is scored with
the same exposure-mapped Poisson rate ratio as single compounds. Each group has its compound count
(`n_compounds`) and the number of members with target counts (`n_with_target`).

Pooling helps because noise and signal average differently. A compound has a handful of counts, but
an A–B disynthon in KinDEL pools up to 190 compounds. MAPK14's target replicates correlate at
r ≈ 0.01–0.08 per compound, and at r ≈ 0.85 once counts are summed per A–B disynthon (the `replicates`
figure in the report). The lines in the cube plot are these disynthons.

Disynthon ratios are **overdispersed**. Synthesis yield and bead stickiness are shared by a whole family,
so a disynthon's 95 % Poisson bound is above 1 for about a quarter of all MAPK14 A–B pairs. A Poisson
interval does not cover family-level variation. This is why the triage thresholds below ask for more than
"bound > 1".

## 6a. Building blocks: a searchable resource (`library/blocks.py`)

KinDEL lists every building block twice: as the free reagent (position C lists a reagent pair) and as the
**reacted fragment**, where pseudo-atoms (francium, astatine) mark the attachment points. The pipeline
publishes one row per block (`blocks/blocks.parquet`, with TSV and Parquet exports and FAIR metadata):

- **Structure.** The canonical fragment with attachment points as `*` (for substructure search), the
  fragment capped with hydrogens, its InChIKey and Murcko scaffold. Capping with H turns a carboxylic acid
  attachment into an aldehyde and so on: descriptors are for comparing blocks, not absolute values.
- **Descriptors** (RDKit): molecular weight, Crippen logP, TPSA, hydrogen-bond donors and acceptors, rotatable
  bonds, heavy atoms, sp3 fraction, rings.
- **Shared fragments.** Blocks whose capped fragment is identical (a methyl and a tert-butyl ester of the
  same acid, or amidine salts) are the same compound part under different DNA codes. In KinDEL,
  21 blocks share a fragment with another block. Their enrichments should agree, which is a free
  consistency check for the data and a hint for library design.
- **Synthesis-yield proxy.** The pre-selection load is one experiment shared by all screens (the pipeline
  checks this). For each block, the mean load of its possible compounds and the share of them that were
  sequenced at all show whether the block is made well. 292 of 1,377 blocks in the sub-library have a mean load
  below half of their position's median (flag `low_load`). A silent block is a yield or truncation suspect,
  not a proven non-binder.

The tables list 377, 1,127 and 190 blocks at positions A, B and C (the paper reports 378, 1,128 and 191).
Position A is limited to the 60 blocks of the sub-library (section 2); B and C are complete.

## 6b. Building-block set enrichment (`analyse/blockset.py`)

The gene-set-enrichment idea, applied to synthons. A building block is a "gene set": the compounds that
contain it. After ranking all compounds of a screen by their Poisson lower bound, a block that drives binding
has its compounds concentrated at the top. For every block the pipeline runs a one-sided Mann-Whitney U test of
its compounds' ranks against all others, computed from rank sums so that 1,377 blocks over 11 M compounds need
one ranking and one group-by. Ties (most compounds score 0) are handled with the standard correction, and the
p-values are adjusted over all blocks (Benjamini-Hochberg). Per block the output holds the AUC (the probability
that one of its compounds outranks a compound without it), the p- and q-values and the number of blocks in the
screen's top 2,000. The test agrees with `scipy.stats.mannwhitneyu` in the unit tests.

Read the AUC, not the q-value. With 11 M compounds nearly every block is "significant" (505 of 1,127 B blocks
in MAPK14), so the q-value stops discriminating. The universe is the compounds seen in at least one column;
compounds absent from all columns are mostly failed syntheses and would otherwise reward blocks that are
made badly.

**A check against known chemistry.** Carbonic anhydrase is inhibited by primary sulfonamides, which bind the
active-site zinc. KinDEL contains exactly three primary-sulfonamide blocks, all at position B. Without being
told, the BCA screen ranks the two para-sulfonamides **first and second** of 1,127 B blocks (AUC 0.996 and
0.990; they are in 1,262 and 596 of the top 2,000 compounds) and the meta isomer 15th. The chance that all
three land in the top 15 at random is about 2 in a million. In the kinase screens the same blocks rank 112th to
799th (AUC 0.45 to 0.59), so the enrichment is specific to the carbonic anhydrase.

The test has limits. Compound sets overlap (every compound is in three sets), so the tests are not independent
and the q-values are indicative. A block-level test finds blocks that dominate a family (a plane of the cube),
not a single disynthon (a line): a 16-compound family inside a 377-compound block does not register, which is
what the disynthon table of section 6 is for.

## 7. Hit triage (`analyse/triage.py`)

The candidates are the `--top_n` compounds (default 2,000) of each **target** screen by Poisson lower bound.
Each candidate gets five flags and a human-readable `reasons` string (`triage/hits.parquet`):

| Flag | Rule | Why |
|---|---|---|
| `replicate_support` | Counts in ≥ 2 of 3 target replicates | One lucky replicate is the classic false positive |
| `disynthon_support` | Best disynthon (AB, AC or BC), **leaving the compound out**, has a rate-ratio bound > 1.5 with ≥ 3 other members with counts | A family supports a compound only if the *other* members are enriched. Without leave-one-out, a single bright compound makes its own disynthons look enriched. |
| `counter_binder` | Bound > 1 in any counter-screen (BCA) | Binds an unrelated protein: promiscuous or an aggregator |
| `matrix_binder` | Bound of matrix-over-load rate ratio > 4 | Enriched on the beads relative to its synthesis yield. The threshold is high because matrix and load are two noisy measurements of the same yield. |
| `shared_target` | Also in the top of another target screen | Two kinases sharing a hinge binder is plausible: **reported, not penalised** |

**Tiers:**
- **1:** replicate and disynthon support.
- **2:** one of the two.
- **3:** neither.
- **flagged:** a counter or bead binder, whatever its support.

Only counter-screen and bead enrichment mark artefacts. Treating "hit in more than one target" as
promiscuity would throw away exactly the selective-within-family chemistry a kinase programme wants.

The thresholds are module constants, not parameters, and are documented in the code. They were chosen on
the KinDEL data, not tuned to the simulator. On real data the results barely depend on them: every DDR1
top hit has disynthon support at bound > 1.5 and still at > 5.

On KinDEL (top 2,000 per kinase, MAPK14 / DDR1):

| | Tier 1 | Tier 2 | Tier 3 | Flagged | of which counter-screen (BCA) | of which bead | Shared target |
|---|---|---|---|---|---|---|---|
| MAPK14 | 1,755 | 75 | 0 | 170 | 168 | 2 | 17 |
| DDR1 | 1,833 | 0 | 0 | 167 | 167 | 0 | 17 |

Almost every top compound by Poisson lower bound already has replicate and family support, because the
bound itself demands evidence. What triage adds is the **counter-screen**: about 8 % of each kinase's top
hits also bind carbonic anhydrase. Of the 17 compounds that are top hits of both kinases, 16 are tier 1.

**What each score lets through.** The same triage, applied to the top 2,000 of every score in turn
(`by_score` in `triage_qc.json`, and a table in the report):

| Ranked by | MAPK14 tier 1 | MAPK14 flagged | DDR1 tier 1 | DDR1 flagged |
|---|---|---|---|---|
| Counts | 1,771 | 150 | 1,811 | 189 |
| Enrichment ratio | 1,672 | 194 | 1,856 | 144 |
| Poisson lower bound | 1,755 | 170 | 1,833 | 167 |
| Faver z (target only) | 250 | **1,744** | 1,499 | 501 |
| EB shrunken ratio | 1,748 | 175 | 1,837 | 163 |

The target-only z-score ranks by target counts alone. In MAPK14, 1,732 of its top 2,000 are **bead
binders**: they are as enriched in the no-protein control as in the target. Any score that divides by the
matrix control removes them before triage sees them.

## 7a. Resynthesis shortlist (`analyse/shortlist.py`)

Hit validation (resynthesis, then a binding assay) is the expensive step of a DEL campaign, and 2,000 triaged
hits are not a decision. Ranking alone also buys redundancy: an enriched disynthon carries all its members to
the top together. The shortlist asks for a diverse plate instead.

1. **Eligible hits**: tier 1 only (replicate and disynthon support, no counter-screen or bead binder).
2. **Cluster** them by structure: Butina clustering on Morgan fingerprints (radius 2, 2048 bits), two hits
   sharing a cluster at Tanimoto similarity of at least 0.65. This is the same family of method KinDEL used to
   choose its validation compounds.
3. **Pick across clusters, best first**: the best hit of every cluster (clusters ordered by their best rank),
   then the second best of every cluster, up to 2 per cluster, and at most 3 per disynthon, until 96 compounds
   (one plate) are chosen. Every number is a parameter.

| 96 compounds | MAPK14 naive top | MAPK14 shortlist | DDR1 naive top | DDR1 shortlist |
|---|---|---|---|---|
| Structural clusters | 17 | **96** | 29 | **96** |
| Disynthons | 5 | **62** | 21 | **58** |
| Mean pairwise similarity | 0.43 | **0.25** | 0.44 | **0.34** |
| Median rank among the top 2,000 | 50 | 441 | 51 | 358 |
| Lowest Poisson lower bound | 11.0 | 2.7 | 35.9 | 14.7 |

Ranking alone buys 96 compounds from 17 families (MAPK14); the shortlist buys 96. The price is stated: the
picks sit deeper in the ranking, but every one of them is enriched by at least 2.7-fold with 95 % confidence,
and the three enriched compounds of a family are each worth a resynthesis less than a new family. Whether a
cluster is a real series or a coincidence of fingerprints is a chemist's judgement; the shortlist makes that
judgement possible on 96 compounds instead of 2,000.

## 8. DESeq2 thinking: what transfers and what breaks

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

## 9. Limitations

- **On-DNA vs off-DNA.** DEL counts measure binding while attached to DNA. Off-DNA KD correlates only loosely,
  as the KinDEL paper also reports.
- **No competition arm.** KinDEL has no selection with a known inhibitor, so allosteric and orthosteric binders
  cannot be separated.
- **Triage thresholds are heuristics.** They are lower bounds of Poisson intervals. Family-level
  overdispersion is handled by asking for more than a bound of 1, not by modelling it. A negative-binomial
  or hierarchical model of synthon effects would turn the thresholds into calibrated probabilities.
- **One counter-screen.** BCA is a single unrelated protein. A compound that binds kinases and nothing
  else on the panel can still be an aggregator.
- **Sub-library.** Results cover about 16 % of BB-As. The pipeline runs on the full library with
  `--n_bb_a 377` given enough memory and time.
- **GPU.** Not used. Every step is vectorised NumPy/SciPy on CPU, and outputs must be byte-stable for snapshot
  tests. A structure-based model (fingerprints → Poisson loss) is a natural next step and could use a GPU.

## References

- Chen B., Danel T., et al. *KinDEL: DNA-Encoded Library Dataset for Kinase Inhibitors.* ICML 2025, arXiv:2410.08938.
- Gerry C. J., Wawer M. J., Clemons P. A., Schreiber S. L. *DNA barcoding a complete matrix of stereoisomeric small molecules.* JACS 141, 10225 (2019).
- Faver J. C., et al. *Quantitative comparison of enrichment from DNA-encoded chemical library selections.* ACS Comb. Sci. 21, 75 (2019).
- Love M. I., Huber W., Anders S. *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.* Genome Biol. 15, 550 (2014).
- Decurtins W., et al. *Automated screening for small organic ligands using DNA-encoded chemical libraries.* Nat. Protoc. 11, 764 (2016).
- Favalli N., Bassi G., Scheuermann J., Neri D. *DNA-encoded chemical libraries – achievements and remaining challenges.* FEBS Lett. 592, 2168 (2018).
