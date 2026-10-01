# STATUS

Handoff note. Written for an assistant session resuming this work cold. Last updated
1 October 2026.

## Where things stand in one paragraph

The manuscript "Avian flight accelerates genome evolution" was rejected by Nature Ecology
and Evolution and then by Nature Communications, six referee reports in total. The project
is being rebuilt for Molecular Biology and Evolution around a single continuous flight
predictor instead of the categorical flight-habit scheme that every referee attacked. This
repository is new and clean. Scripts 01 to 13 are carried over unchanged and work.
Scripts 14 and 15 are written and tested. Scripts 16 to 24 are specified but not written.

## Paths

| | |
|---|---|
| This repo | `~/Dropbox/Research/auk/avian_flight_genomics` |
| Old repo, do not write to it | `~/Dropbox/Research/auk/avi_family_flight_genomics` |
| Full plan | `~/Dropbox/Research/auk/PLAN_next_submission.md` |
| Submission-ready table | `~/Dropbox/Research/auk/Submission_MBE/Table_S1.xlsx` (identical to `data/Table_S1.xlsx`) |
| Referee reports | `~/Dropbox/Research/auk/Decision on Nature Ecology...txt` and `Submission_NatComms/Decision on manuscript NCOMMS...txt` |
| AVONET | `~/Desktop/Research/predicts_colour/avonet/AVONET Supplementary dataset 1.xlsx`, sheet `AVONET1_BirdLife` |
| Livezey 1988, Poole 1938 | `~/Dropbox/Research/auk/Livezey.pdf`, `Poole.pdf` |
| Codon alignments and gene trees | `~/Desktop/Research/auk/postshrink_alignments` (not in the repo) |
| HyphAeon notes | `HYPHAEON.md` in this repo |

## Decisions already taken. Do not reopen these.

1. **84 taxa, 42 sister pairs, at every stage.** The user ruled out extending to the 363
   B10K taxa twice. Do not propose it again, even though AVONET would cover 358 of them.
   The reason is that flight ability and primary flight style are hand-curated for 84 only.
2. **Primary predictor is `ln(mass / wing length^2)`**, a wing-loading index, replacing the
   four-level flight-habit categories.
3. **Natural log throughout**, matching `masdif`.
4. **HyphAeon `phenotype` for the per-gene work**, with RERconverge as the peer-reviewed
   parallel. RELAX and BUSTED-PH were considered and rejected as not fitting a
   rate-association framework.
5. **GSEA on the continuous per-gene statistic**, never a top-percentage cut. The previous
   submission died partly on an over-representation analysis at a 10% cut with q < 0.5.
6. **MR-PMM replaces the Bayesian path analysis.** `brms` is installed; `MCMCglmm` is not.
7. **Flight categories kept as published.** AVONET disagreed with the user's migration
   coding on only 8 of 83 taxa and the user's coding won 6 of those, so the categories
   stand.

## Numbers already established. Reuse rather than recompute.

**Why the index.** Per taxon, variance explained:

| variable | flight ability | flight style | collinear with migration |
|---|---|---|---|
| HWI raw | 0.238 | 0.424 | yes |
| ln(HWI/mass), old variable | 0.334 | 0.514 | no, but 78% body mass |
| ln(wing / mass^(1/3)) | 0.490 | 0.682 | no |
| **ln(mass / wing length^2)** | 0.453 | **0.683** | no |

At contrast level all four land between R² 0.39 and 0.46, so the index's advantage is at
the taxon level, not in contrasts. Only ΔHWI is confounded with migration (p = 0.020);
the index is not (p = 0.300).

**Validation.** Published wing areas exist for 8 of the 84 taxa. Across them the index
tracks measured wing loading at R² = 0.93 (0.92 if each source's own specimen mass is used on both sides), ratio 0.67 to 2.14, spanning 17 g to 4664 g.
Wing area scales as wing length^1.93 (95% CI 1.51 to 2.35), which is why the exponent is 2.

**Design arithmetic after relabelling.** minimal 15, to migratory 15, toward flight loss 7,
toward limited flight 5. Five of the seven flight-loss pairs also carry a
sedentary-to-migratory change, so only pairs 1 and 13 isolate a single axis.

**Prior results worth keeping in view.** The only significant genome-wide signal is on the
migration axis, dS and dN in sedentary-to-migratory pairs at permutation p = 0.028 to
0.046. Flight loss was null in all three metrics, p > 0.4. An analysis built only on the
flight index may well return null, which is why script 16 fits migration alongside it.

**The great auk.** Wing length 160 ± 9 mm (n = 23 skins, Livezey 1988 Table 1), wing area
230 cm², body mass 3561 g (Montgomerie and Birkhead 2025, replacing Smith 2016's PGLS
estimate of 4750 g), wing loading 15.5 g/cm², which is 9.5 times *Alca torda*'s 1.63. Its
HWI of 53.84 is razorbill-like, so it kept the pointed wing *shape* while halving area and
multiplying mass. On the index it is second lowest of all 84 taxa, above only
*Aptenodytes forsteri*.

## Scripts

Done and tested: `01` to `13` (carried over), `14` to `22` except `23`.
To write: `24_figures`.

**RERconverge is not needed and is not used.** HyphAeon's gene p-value already uses the
RERconverge null (Brownian-motion permulations), so installing the package would supply no
new null. What it would supply is a different *statistic*, and that is worth having, since
HyphAeon is mammal-trained with avian transfer benchmarked only on mitochondrial data.
Script 20 therefore computes the comparator directly from `dnds_branches.Rdata`: relative
evolutionary rates recovered by double-centring the log per-locus per-taxon rates, which
is the residual RERconverge models, then correlated with each trait. 13,625 genes in
minutes, from PAML rates we generated ourselves. No external package, no treesObj build.

`clusterProfiler` 4.18.4 and `fgsea` 1.36.2 are now installed. Script 21 builds gene sets
from `goterms.csv` directly (3,025 GO terms, 389 within the 10 to 500 size bounds), so no
OrgDb is needed and the annotation is identical to the previous submission.

**`23_mrpmm` does not converge as specified.** A short run gave max Rhat 1.99, minimum ESS
ratio 0.01 and 40 divergent transitions, so the correlations it prints must not be
interpreted. Separating phylogenetic from residual covariance across several responses is
weakly identified with 81 taxa. The script now defaults to four responses rather than six,
coerces the tree to ultrametric (it fails `is.ultrametric` by rounding, which left the
covariance badly conditioned), uses a tighter prior on the phylogenetic sd, and refuses to
let a bad fit pass silently by warning when Rhat exceeds 1.05. A production run is
pending. If it still fails, the options in order are three responses, or fixing the
phylogenetic signal rather than estimating it per response.

`18` and `19` are tested: 14,124 of 14,957 genes retained (median 77 phenotyped taxa), and
the runner completed four genes with no failures at 0.9 to 3.2 s each. At that rate one
phenotype is about 9 hours and all three about 26. Each is specified in `PLAN_next_submission.md`.

`20_hyphaeon_collect` must also correlate the per-gene HyphAeon statistic against
`aitchison_variance_results.csv`. That file is a per-gene measure of base-composition
heterogeneity and exists to answer the referee asking whether genes are being selected
for downstream analysis by their base composition. A null correlation is the result
being sought. This is a per-gene check and is unrelated to the per-taxon GC covariate
that script 23 uses.

Script 14 writes `output/traits_per_taxon.csv` (84 rows) and `output/traits_contrasts.csv`
(42 rows). Everything downstream reads those two files, not the spreadsheet.

## Results so far

**16, rate contrasts.** The flight index explains nothing that survives BH correction
across the six responses (all p_adj = 0.123; raw p for dS is 0.049 to 0.055). Migration
does: omega is lower in the migratory sister, p_adj = 0.019 on both response types. This
is the pattern anticipated before the script was written, and it is why migration is in
the model.

**17, power.** With 42 pairs, 80% power needs a correlation of 0.41 for a single gene at
alpha = 0.05, and 0.55 for a gene inside a BH-controlled screen of 14,785. That is a hard
bound on what scripts 18 to 20 can return, and a null there must be read against it.

**22, age.** The age-by-trait interaction is not significant for any rate type
(p = 0.68, 0.99, 0.18). The age effect sits in the minimal-change pairs (dN slope
p = 0.023, dS p = 0.077) and is absent from the pairs that changed (p = 0.99, 0.70). The
latency simulation cannot generate the observed age correlation at any latency for dN or
omega. **The previous paper's reading of its Figure 3, that older pairs show larger
trait-driven differences, does not survive.** The age relationship is elapsed time, not
trait effect.

## Traps

- **`sp1` in `dnds_sisdifs_vol.Rdata` is NOT always the lesser-flight sister.** In 7 of
  the 42 pairs it is the greater-flight one: ANTCAR-CHOACU, AQUCHR-HALALB, STROCC-TYTALB,
  TROMEL-APAVIT, BUCABY-BUCRHI, MENNOV-ATRCLA, TAEGUT-LONSTR. These are the pairs AviList
  reordered. `moldiff` is sp1/sp2, so orienting by position silently flips the sign for a
  sixth of the data. Script 16 orients against the trait table and prints the affected
  pairs. Any new script touching `moldiff` must do the same.
- **Species names differ between sources.** Table S1 uses AviList 2025, the tree and
  `taxonomy_codes.csv` use B10K names. Join through the six-letter taxon code, never by
  name. Known renames: *Laterallus rogersi* = `Atlantisia_rogersi`, *Anarhynchus
  alexandrinus* = `Charadrius_alexandrinus`, *Nannopterum harrisi* = `Phalacrocorax_harrisi`,
  *Crinifer concolor* = `Corythaixoides_concolor`, *Dryobates pubescens* =
  `Picoides_pubescens`. `PINIMP` is absent from `taxonomy_codes.csv` and is patched in by
  script 14.
- **AVONET synonyms**, needed only if querying AVONET directly: *Crinifer concolor* ->
  *Corythaixoides concolor*, *Anarhynchus alexandrinus* -> *Charadrius alexandrinus*,
  *Alaudala cheleensis* -> *Alaudala rufescens* (a genuine split, not a pure synonym).
  *Pinguinus impennis* is absent from AVONET entirely.
- **The old `data/Table_S1.xlsx` carried a body mass of 3560 g for *Laterallus rogersi*,
  about ninety times too high**, copied from the *Pinguinus* cell. This repo's copy is
  corrected. If you ever see 3560 against *Laterallus*, you are reading a stale file.
- **`masdif` in `genome_flight_age_mass_data.Rdata` for ATLROG-HELFUL does not reconcile
  with any version of the trait table.** Stored value 1.5126; the current clade masses give
  0.8924. This pair is also the most influential observation in the old path analysis.
  Unresolved. Track down before trusting any model that uses `masdif` for that pair.
- **Excel shared formulas.** When editing `Table_S1.xlsx` by zip surgery, a regex matching
  `<f>...</f>` catches only the plain formulas. 82 of the 168 cells use
  `<f t="shared" si="n"/>` and reference a master cell. Patch by column and row range, not
  by formula text, or you will silently leave half the sheet on the old scale.
- **`migration_pair` is 21 pairs, not 15.** It follows the Migration column of Table S1
  (`sed.->migr.`), whereas the 15 is the count of the "to migratory" flight-habit category.
  They differ because several flight-loss pairs also carry a migration change. This is
  deliberate and is the confound worth reporting.
- `02_align_basic_iq` sources `trimCols.R` and `trimCodons.R`, which were missing from the
  old repository and are now in `functions/`.

## Not installed

`clusterProfiler` (Bioconductor, needed by 21), `RERconverge` (GitHub, needed by 20),
`ggtree`, `MCMCglmm`. Available and used: `readxl`, `ape`, `phangorn`, `ggplot2`, `dplyr`,
`brms`, `phytools`, `patchwork`. `node` is not installed, so the dataviz palette validator
cannot be run here.

## HyphAeon: resolved, see HYPHAEON.md

Installed and run to completion on a real gene. All four questions answered.

1. **It is RERconverge-style.** The `--permulations` flag is documented as the
   "RERconverge null model", so the gene-level p-value is an empirical permulation p-value
   under Brownian motion, not a diversifying-selection test. The plan stands unchanged.
2. **Inputs match our files unchanged.** 100% taxon concordance, 76 shared taxa on the test
   gene. Alignments and per-gene trees are at `~/Desktop/Research/auk/postshrink_alignments`
   (about 14,957 genes). Taxon labels are the six-letter codes.
3. **1.5 s of compute per gene**, but 5.3 s per process, so the runner must load the model
   once and loop rather than calling the CLI per gene.
4. **Gene trees**, which sit beside each alignment and are accepted directly.

**Gene-level statistics are in the JSON only**; the CSV is site-level. The collector must
parse JSON.

**Environment, Intel Mac.** torch caps at 2.2.2, which predates NumPy 2, so the whole stack
must be pinned to the NumPy 1 era *before* installing hyphaeon. The exact recipe is in
HYPHAEON.md. Also recorded there: `exceptiongroup` segfaults inside its own error
formatter, hiding every real error until the hook is reset.

## Open questions, blocking

1. **What are the input and output formats of `hyphaeon phenotype`?** Needed to write the
   writer in script 18 and the parser in script 20. HyphAeon 0.2.0 is on PyPI, so a `pip
   install hyphaeon` and `hyphaeon phenotype --help` settles this and the two below.
2. **Is the HyphAeon per-gene statistic a rate association, as in RERconverge, or a
   diversifying-selection test, as in MEME?** The preprint names PARS as a pillar and gives
   no methods for it. This matters: the user's reading of a genome-wide excess as relaxed
   selection holds only under the first interpretation, because relaxation moves omega
   toward 1 from below where a diversifying test has little power.
3. Runtime and batching for roughly 44,000 gene analyses (3 phenotypes x ~14,800 genes).
4. Gene trees or the dated species tree per gene analysis.

**Great auk damage control: partly already done.** `Submission_NatComms/figS1` IS the
mapDamage output. It shows base composition at read ends, misincorporation frequencies by
position, read-length distributions (mode near 45 bp) and cumulative C-to-T and G-to-A
curves. Damage is low, roughly 0.02 to 0.04 at the termini. It was added for the Nature
Communications submission after NEE Reviewer 1 asked for it.

What is still missing is contamination inference and, more importantly, a quantification
of how the residual deamination shifts dN, dS and omega for *Pinguinus* specifically.
That sensitivity check should run before *Pinguinus* is trusted in any codon-level result,
but it no longer blocks starting script 19.

## Next action

Write `23_mrpmm.Rscript`, which is unblocked and uses `brms`. Scripts 18 to 21 wait on the
two HyphAeon questions above, and 19 additionally on the great auk damage control.
