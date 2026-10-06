# STATUS

Handoff note. Written for an assistant session resuming this work cold. Last updated
6 October 2026.

## Where things stand in one paragraph

The manuscript "Avian flight accelerates genome evolution" was rejected by Nature Ecology
and Evolution and then by Nature Communications, six referee reports in total. The project
is being rebuilt for Molecular Biology and Evolution around a single continuous flight
predictor instead of the categorical scheme every referee attacked. Scripts 01 to 26 are
written and run. The per-gene screens are finished, the MR-PMM has been fitted in four
specifications and compared, and the figures are drawn. The open work is a per-site PELICAN
screen, which cannot run on this Mac and is packaged for GenomeDK in `hpc/`.

## Paths

| | |
|---|---|
| This repo | `~/Dropbox/Research/auk/avian_flight_genomics` |
| Old repo, do not write to it | `~/Dropbox/Research/auk/avi_family_flight_genomics` |
| Full plan | `~/Dropbox/Research/auk/PLAN_next_submission.md` |
| Submission table | `~/Dropbox/Research/auk/Submission_MBE/Table_S1.xlsx` (same as `data/Table_S1.xlsx`) |
| Referee reports | `~/Dropbox/Research/auk/Decision on Nature Ecology...txt`, `Submission_NatComms/Decision on manuscript NCOMMS...txt` |
| AVONET | `~/Desktop/Research/predicts_colour/avonet/AVONET Supplementary dataset 1.xlsx`, sheet `AVONET1_BirdLife` |
| Livezey 1988, Poole 1938 | `~/Dropbox/Research/auk/Livezey.pdf`, `Poole.pdf` |
| Codon alignments and gene trees | `~/Desktop/Research/auk/postshrink_alignments` (not in the repo) |
| HyphAeon venv (Intel Mac pinned build) | scratchpad `hvenv3`; recipe in `HYPHAEON.md` |

## Decisions already taken. Do not reopen these.

1. **84 taxa, 42 sister pairs, at every stage.** Ruled out extending to the 363 B10K taxa
   twice. Do not propose it again, even though AVONET would cover 358 of them. Flight
   ability and primary flight style are hand-curated for 84 only.
2. **Primary predictor is `ln(mass / wing length^2)`.**
3. **Natural log throughout**, matching `masdif`.
4. **HyphAeon `phenotype` for the per-gene work.** Confirmed RERconverge-style.
5. **RERconverge the package is NOT used.** It would add no new null. The comparator is
   computed in script 20 from `dnds_branches.Rdata` instead.
6. **GSEA on the continuous statistic**, never a top-percentage cut.
7. **MR-PMM replaces the Bayesian path analysis.** `brms`, not `MCMCglmm`.
8. **Flight categories kept as published.** AVONET disagreed on only 8 of 83 taxa and our
   coding won 6 of those.
9. **Residual correlation is OFF in script 23.** It was retained through two fits, both
   of which failed to mix. The failing version is kept as evidence under `*_rescor` and
   is reproducible with `MRPMM_RESCOR=TRUE`.

## Numbers already established. Reuse rather than recompute.

**Why the index.** Per taxon, variance explained:

| variable | flight ability | flight style | collinear with migration |
|---|---|---|---|
| HWI raw | 0.238 | 0.424 | yes |
| ln(HWI/mass), old variable | 0.334 | 0.514 | no, but 78% body mass |
| ln(wing / mass^(1/3)) | 0.490 | 0.682 | no |
| **ln(mass / wing length^2)** | 0.453 | **0.683** | no |

At contrast level all four land between R² 0.39 and 0.46, so the advantage is at the taxon
level. Only ΔHWI is confounded with migration (p = 0.020); the index is not (p = 0.300).

**Validation.** Wing areas exist for 8 of 84 taxa. The index tracks measured loading at
R² = 0.93, ratio 0.67 to 2.14, across 17 g to 4664 g. Wing area scales as wing
length^1.93 (95% CI 1.51 to 2.35), which is why the exponent is 2.

**Design arithmetic.** minimal 15, to migratory 15, toward flight loss 7, toward limited
flight 5. Five of the seven flight-loss pairs also carry a migration change, so only pairs
1 and 13 isolate a single axis.

**The great auk.** Wing length 160 ± 9 mm (n = 23 skins, Livezey 1988 Table 1), wing area
230 cm², mass 3561 g (Montgomerie and Birkhead 2025, replacing Smith 2016's 4750 g), wing
loading 15.5 g/cm², 9.5 times *Alca torda*'s 1.63. Its HWI of 53.84 is razorbill-like, so
it kept the pointed wing shape while halving area and multiplying mass. Second lowest of
all 84 on the index, above only *Aptenodytes forsteri*.

## Results so far

**16, rate contrasts** (corrected after the factor-indexing bug; the earlier numbers in
this file were wrong). The flight index explains nothing surviving BH correction, minimum
p_adj 0.168. Migration does, and the corrected result **reproduces the original paper's
published finding**: dN and dS are elevated in the migratory sister (p_adj 0.0078 and
0.015, positive estimates) while omega is unaffected (0.24, 0.18). The buggy version had
claimed the opposite, a significant omega effect and null dN and dS.

**20b, contrast screen** (also corrected). One gene passes q < 0.05, GALGAL_R08885 on
dN against the flight index (n = 36, r = -0.695, q = 0.017). The seven genes reported
before the fix were artefacts.

The informative result is genome-wide rather than per gene. Against a 5% expectation:

| rate | phenotype | genes below p = 0.05 | median p |
|---|---|---|---|
| dS | migration | **11.7%** | 0.377 |
| dN | migration | **9.5%** | 0.412 |
| dS, dN, omega | loading | 4.2 to 4.6% | 0.52 to 0.54 |
| dS, dN, omega | flightloss | 2.7 to 3.1% | 0.55 to 0.57 |

So migration produces a strong, genome-wide excess of small p-values on dS and dN, matching
script 16 and the original paper, but no single gene survives correction because the binary
test's resolution floor (minimum p 1.3e-05 against a 4.7e-06 requirement) blocks it. The
flight index sits at the nominal rate, that is, no signal.

**17, power.** 80% power needs rho = 0.41 for a single gene, **0.55 inside a BH-controlled
screen of 14,785**. A null in the screen must be read against that bound.

**22, age.** The age-by-trait interaction is null for every rate type (p = 0.68, 0.99,
0.18). The age effect sits in the minimal-change pairs (dN slope p = 0.023) and is absent
from pairs that changed (p = 0.99). The latency simulation cannot generate the observed
age correlation at any latency for dN or omega. **The previous paper's reading of its
Figure 3 does not survive.** The age relationship is elapsed time, not trait effect.

**Variable audit for script 23** (per taxon, n = 83):

| response | skew | kurtosis | lambda |
|---|---|---|---|
| ln dS | −0.32 | −0.37 | 1.03 |
| ln omega | 1.27 | 3.04 | 0.00 |
| **loading index** | **1.67** | **7.08** | 0.96 |
| ln genome length | −0.62 | 0.94 | 0.60 |
| ln gene count | −0.09 | −0.14 | 0.31 |
| GC content | 1.12 | 2.62 | 0.91 |
| ln mean intron length | −0.83 | −0.26 | 0.00 |

The flight index is the one badly skewed variable, driven by the wing-propelled divers and
the ratites. Real biology, but a Gaussian response is an approximation there. Max pairwise
|r| among responses is 0.574.

## Scripts

Done and tested: `01` to `22`, including the new `20b_contrast_screen`, which is now
the primary per-gene analysis. To write: `24_figures`.

Running now: `19_hyphaeon_run.py` for all three phenotypes in parallel, 14,124 genes each,
1000 permulations, `--resume`. At about 69% after roughly 19 hours; the rate settled to
about 7 s per gene rather than 5.07 because three processes contend on 16 cores. Roughly
8 more hours. Then `20`, `21`.

**The flightloss phenotype has 92 legitimate failures**, all `Insufficient foreground taxa
(1)`. The binary flight-loss phenotype has only 7 positives across 84 taxa, so genes with
sparse taxon sampling contain just one and cannot support a foreground/background test.
That is the n = 7 limit of the design showing up gene by gene, not a bug. loading and
migration have zero failures.

`23_mrpmm` is on its **third** specification, and the first two are documented failures.

| version | responses | rescor | result |
|---|---|---|---|
| 1 | ln dN, ln dS, ln len, flight index | TRUE | Rhat 1.99, 40 divergences, an hour without finishing |
| 2 | dS, omega + 5 more (7 total) | TRUE | Rhat 1.45, ESS ratio 0.007, 25 divergences, low BFMI |
| 3 | same 7 | **FALSE** | running |

Version 1 failed because ln dN and ln dS correlate at r = 0.970 and lambda was at or above
1 for three of four responses. Replacing them with **dS and omega** cut collinearity to
0.505 and gave omega lambda 0.00, which helped but was not enough.

Version 2's low BFMI is the diagnostic that matters: the sampler cannot move through the
region where phylogenetic and residual variance trade off. Lambda is 1.03 for ln dS, 0.96
for the flight index and 0.91 for GC, so three of seven responses have essentially no
non-phylogenetic variance, and 21 residual correlations were being estimated from variances
pinned at the floor.

**Version 2 is kept as evidence**, by decision, under the `*_rescor` suffix:
`output/mrpmm_fit_rescor.Rdata`, `mrpmm_correlations_rescor.csv`,
`mrpmm_diagnostics_rescor.csv`, `mrpmm_rescor_run.log`,
`figures/Fig3_mrpmm_rescor.pdf`. Reproduce with `MRPMM_RESCOR=TRUE`.

Version 3 drops `set_rescor`, removing 21 parameters the data cannot supply. The cost,
which belongs in the paper rather than buried in a comment, is that the split between
phylogenetic and residual covariance becomes partly assumed rather than estimated. If it
still fails, cut the response set to three or four.

## Timings, all measured

| script | time |
|---|---|
| 14, 15, 16, 22 | 1 to 6 s each |
| 17 power | 28 s |
| 18 prep | 33 s |
| **19 screen** | **5.07 s per gene, ~20 h per phenotype, ~20 h wall for three in parallel** |
| 20 collect | ~15 to 30 min at full scale |
| 21 GSEA | ~2 to 5 min at full scale |
| 23 MR-PMM | unresolved, over an hour on the failing specification |

Permulations are free: the same gene at 0, 100 and 1000 took 12.7, 11.3 and 11.5 s.


## Per-gene screening: what is appropriate, and why the alternatives are not

Settled empirically. Full reasoning in `HYPHAEON.md`; the short version follows.

**The primary screen is `20b_contrast_screen`**, correlating each gene's sister-pair rate
contrast against the contrast in the continuous flight index across 42 pairs. Sister pairs
are independent contrasts by construction, so no phylogenetic model, no lambda and no
permulation null enter, and there is nothing to misspecify. Median p is 0.496 to 0.497,
and it finds 6 genes on dS and 1 on omega at q < 0.05. Script 17 computed power for
exactly this test.

**Only two of nine rate-by-phenotype screens are usable**, and both are the continuous
flight index. The rest fail a *resolution floor*: a discrete test cannot report a p below
its own minimum, and if that floor sits above the Benjamini-Hochberg threshold for the top
gene the screen cannot make a discovery whatever the biology.

| screen | distinct p | min p | BH needs | usable |
|---|---|---|---|---|
| dS x loading | 10,698 | 3.3e-07 | 4.7e-06 | yes |
| omega x loading | 6,908 | 3.3e-06 | 7.2e-06 | yes |
| dS x flightloss | 2,673 | 8.4e-06 | 4.8e-06 | no |
| dS x migration | 2,982 | 6.2e-04 | 4.7e-06 | no, 130-fold short |
| dN x loading | 6,719 | 6.9e-05 | 7.4e-06 | no |
| four others | | | | no |

The binary phenotypes are resolution-limited because 42 pairs split 7/35 or 21/21 give
only a few thousand distinguishable outcomes. This is a property of the design, not of the
test: switching from a t approximation to an exact Wilcoxon made it slightly worse, not
better. dN fails despite being continuous because it has usable data for only 6,719 genes.
The remedy would be a graded phenotype, not a different test, and we have chosen not to
pursue it.

**HyphAeon is not used for inference.** Three reasons, each measured:

1. Its permulation null assumes pure Brownian motion. It is exactly calibrated at
   lambda 0.994 (median p 0.500, 5.0% below 0.05, KS p 0.60) and becomes steadily more
   conservative as the trait departs from Brownian: 0.789 at lambda 0.85, 0.957 at
   lambda 0. Our flight index sits at 0.85, which predicts the whole-screen median of
   0.786 and the zero discoveries. A lambda-matched null would patch this, but trait
   conservation varies across the tree so no single lambda is right anywhere.
2. None of its gene summaries is null-centred (fraction positive from 0.06 to 0.998), so
   none can rank genes for GSEA. An early run on one returned 28 sets at q < 0.05 with
   every NES positive, which is an uncentred list rather than a finding.
3. Rank concordance with relative rates from our own PAML branch lengths is 0.009 to
   0.054, essentially zero, and the comparator found 9 genes where HyphAeon found none.

The full HyphAeon screen is retained in `output/hyphaeon_out` and reported as an
exploratory layer plus a methods observation worth sending to the authors.

**GSEA** ranks on the contrast statistic and reports nothing at q < 0.05. The comparator
ranking yields 6 sets, three of which (GPCR activity, GPCR signalling, olfactory receptor
activity) are a single overlapping, paralogy-prone family and should be treated with care.

## Traps

- **`comparison` in `dnds_sisdifs_vol.Rdata` is a FACTOR, and indexing a named vector with
  a factor uses its integer level codes, not its labels.** This silently paired every
  gene's rate contrast with the wrong pair's trait value in both `16_rate_contrasts` and
  `20b_contrast_screen`, with no error and plausible output. It inverted the script 16
  conclusions and manufactured seven spurious genes in the screen. Always
  `as.character()` before indexing. Both scripts are fixed and carry a comment saying why
  the conversion is load-bearing.
  This is the second silent mis-pairing bug here, after the `sp1` orientation trap. Both
  came from joining pair-level data by position or implicit conversion rather than by an
  explicit key. Any new script touching pair-level data should assert the join.

- **`sp1` in `dnds_sisdifs_vol.Rdata` is NOT always the lesser-flight sister.** In 7 of 42
  pairs it is the greater-flight one: ANTCAR-CHOACU, AQUCHR-HALALB, STROCC-TYTALB,
  TROMEL-APAVIT, BUCABY-BUCRHI, MENNOV-ATRCLA, TAEGUT-LONSTR. `moldiff` is sp1/sp2, so
  orienting by position silently flips the sign for a sixth of the data. Script 16 orients
  against the trait table and prints the affected pairs. Any new script touching `moldiff`
  must do the same.
- **`Feng_Table_1.csv` carries four trinomials.** `Aquila chrysaetos canadensis` and
  `Strix occidentalis caurina` are among them, and matching on the full name silently drops
  those taxa. Match on the first two words, and apply the `Nannopterum harrisi` to
  `Phalacrocorax_harrisi` synonym. That gives 83 of 84.
- **Species names differ between sources.** Table S1 uses AviList 2025, the tree and
  `taxonomy_codes.csv` use B10K names. Join through the six-letter taxon code. Renames:
  *Laterallus rogersi* = `Atlantisia_rogersi`, *Anarhynchus alexandrinus* =
  `Charadrius_alexandrinus`, *Nannopterum harrisi* = `Phalacrocorax_harrisi`,
  *Crinifer concolor* = `Corythaixoides_concolor`, *Dryobates pubescens* =
  `Picoides_pubescens`. `PINIMP` is absent from `taxonomy_codes.csv` and patched in by
  script 14.
- **AVONET synonyms**, if querying AVONET directly: *Crinifer concolor* ->
  *Corythaixoides concolor*, *Anarhynchus alexandrinus* -> *Charadrius alexandrinus*,
  *Alaudala cheleensis* -> *Alaudala rufescens* (a genuine split). *Pinguinus* is absent.
- **The great auk has no comparable genome length.** `genome_size_data.csv` lists 19 Mb for
  PINIMP, which is reads mapped to the razorbill, not an assembly. Exclude it from genome
  length and gene count, keep it everywhere else.
- **The old `data/Table_S1.xlsx` had *Laterallus rogersi* at 3560 g**, ninety times too
  high, copied from the *Pinguinus* cell. This repo's copy is corrected. Seeing 3560 against
  *Laterallus* means a stale file.
- **`masdif` for ATLROG-HELFUL does not reconcile** with any version of the trait table.
  Stored 1.5126; current clade masses give 0.8924. Also the most influential observation in
  the old path analysis. Unresolved.
- **Excel shared formulas.** Editing `Table_S1.xlsx` by zip surgery with a regex matching
  `<f>...</f>` catches only the plain formulas. 82 of 168 cells use `<f t="shared" si="n"/>`
  and reference a master. Patch by column and row range, not by formula text.
- **`migration_pair` is 21 pairs, not 15.** It follows the Migration column
  (`sed.->migr.`), whereas 15 is the "to migratory" flight-habit category. Deliberate, and
  the confound worth reporting.
- **`exceptiongroup` segfaults inside its own error formatter** in the HyphAeon venv, hiding
  every real error. Reset `sys.excepthook` after importing it. Script 19 already does.
- `02_align_basic_iq` sources `trimCols.R` and `trimCodons.R`, now in `functions/`.

## Open questions

1. **Does the revised script 23 converge?** If not, see the ordered options above.
2. **The flight index is skewed** (1.67, kurtosis 7.08). Leave it and check the posterior
   predictive, or use a heavier-tailed response for that variable.
3. **Where to run the screen long term.** It works here with the pinned stack, but a Linux
   or Apple Silicon machine would need no pinning and would be faster.
4. **Great auk contamination inference and a damage sensitivity check.** The mapDamage
   plots already exist as `Submission_NatComms/figS1`, showing low damage (0.02 to 0.04 at
   read termini). What is missing is contamination inference and a quantification of how
   residual deamination shifts dN, dS and omega for *Pinguinus*.

## Carried over, not yet scheduled

GC-biased gene conversion and GC3, effective population size from Stiller et al. 2024,
clutch size and generation length, latitude, and pipeline-wide multiple-testing accounting.

## Next action

Wait for the screen, then run `20_hyphaeon_collect` and `21_gsea`. Check whether the
revised `23_mrpmm` converged. Write `24_figures` last, once the figures it assembles exist.


## Added 6 October 2026

### MR-PMM, settled

Four fits, 9.2 hours each, in `output/mrpmm_fit_norescor*.Rdata`. `scripts/25_mrpmm_compare.Rscript`
reports across them and caches what it needs into `output/mrpmm_harvest_*.csv`, so the
76 Mb fit objects can be deleted without breaking the comparison.

Four correlations are credible in all four and agree to within 0.05. Neither flight axis
correlates with anything. Residual correlation stays off, and trimming outliers makes the
case for it worse rather than better, which is measured in the header of script 25.

**Read convergence by parameter block, not by headline Rhat.** Every failing parameter is
a residual standard deviation or sits on ln dS. The 21 reported correlations are clean in
all four fits.

### PELICAN, packaged but not run

The one resolution nothing here has examined is the individual site. A few sites under
strong selection move neither genome-wide nor gene-level omega, so our null on omega says
nothing about them.

**It cannot be installed on this Mac.** Homebrew refuses to build opam because the Command
Line Tools predate Xcode 26.3, and Docker Desktop 20.10.6 starts its backend but never
opens a socket. Both need `sudo` and a system update. Do not spend another hour on this.

`hpc/` is a complete GenomeDK kit: Apptainer pull, two smoke tests, a Slurm array and a
collector. `hpc/inputs/` is tracked, so a clone carries everything except the alignments,
which `hpc/make_transfer.sh` packs into a 113 Mb tarball.

**The unresolved question is in `hpc/00_setup.sh`.** 7,087 of the 9,890 eligible genes have
a taxon set no other gene shares, so if PELICAN insists that tree and alignment agree,
batching is impossible and every gene needs its own pruned tree. The smoke test decides it
and writes `pelican_mode.txt`. Do not guess, and do not skip it.

PELICAN 1.0.8 has `scan continuous`, so the loading index needs no binarisation. Run order:
flight (discrete, 15 tips), then loading (continuous, costlier), then migration (29 tips,
best powered and worst motivated, since our result there is a mutation-rate signature).

### Trap recorded

**Do not edit a script while Rscript is executing it.** R parses lazily, so an edit lands
mid-parse and kills the run. This cost the post-fit reporting of the trimmed fit, though
the fit object had already been saved.
