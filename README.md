# avian_flight_genomics

Scripts and data for the analysis of avian flight habits and genome evolution, prepared
for submission to Molecular Biology and Evolution.

This repository replaces `avi_family_flight_genomics`. It carries forward the alignment,
dating and sister-pair pipeline unchanged, and replaces the categorical flight predictor,
the gene-set enrichment approach and the Bayesian path analysis. Nothing superseded has
been copied across, so a release from here contains no mixture of old and new material.

`STATUS.md` is the working handoff note: decisions taken, numbers established, traps and
the next action. `HYPHAEON.md` documents the per-gene method, its interface and the install
recipe. Read both before resuming work.

## Design

42 sister pairs, 84 taxa, roughly 14,800 protein-coding genes. All analyses are restricted
to these 84 taxa, which are the species with curated flight traits. No analysis uses the
wider 363-taxon B10K set.

The primary flight predictor is **ln(mass / wing length²)**, a wing-loading index,
replacing the four-level flight-habit categories. Published wing areas exist for only 8 of
the 84 taxa, so area cannot be a study-wide variable. Wing length is available for all 84,
from AVONET for 83 and Livezey (1988, Table 1) for *Pinguinus impennis*. Wing area scales
as wing length to the power 1.9 to 2.4 where both are known, so the index has the
dimensions of wing loading and tracks measured loading at R² = 0.93 across 17 g to 4664 g
and six flight styles. The index uses the AVONET mass while the published loading uses each
source's own specimen mass, so part of the scatter in that comparison is a difference in
mass source rather than in the length² approximation.

Two further axes enter as binary phenotypes: gain of migration, and flight loss. The
wing-loading index is statistically orthogonal to migration, so both can be fitted without
collinearity.

All logarithms are natural, matching `masdif` in `genome_flight_age_mass_data.Rdata`.

## Pipeline

Scripts run from `scripts/`, with paths `../data/`, `../functions/`, `../figures/` and
`../output/`. Script 14 writes `output/traits_per_taxon.csv` (84 rows) and
`output/traits_contrasts.csv` (42 rows); every downstream script reads those rather than
the spreadsheet.

### Carried forward unchanged

| | |
|---|---|
| `01_cut_exons` | Extract exons from B10K c123 data |
| `02_align_basic_iq` | MACSE alignment and IQ-TREE |
| `03_clean_fast_tips` | Remove outlier long branches, top 1% |
| `04_dnds_runs` | PAML and EraBLE for dN and dS trees |
| `05_dnds_collect` | Sister-pair rate ratios |
| `06_blast_loci` | BLAST locus identification |
| `07_id_focal_loci` | Locus bookkeeping. Its former role feeding threshold-based enrichment is replaced by GSEA in `21_gsea` |
| `08_prep_dating` | Select 313 clock-like loci, concatenate |
| `09_binomtests` | Per-gene binomial and KS tests |
| `10_aitchison_variance` | Per-gene base-composition heterogeneity |
| `11_base_composition_vs_rates` | Aitchison variance against binomial p-values |
| `12_ks_permutation_envelopes` | Permutation null envelopes for the KS tests |
| `13_genome_size` | Genome metrics and assembly-quality correction |

### New for this submission

| | status | |
|---|---|---|
| `14_traits_assemble` | run, 2 s | Per-taxon and per-pair trait tables, the wing-loading index, its contrasts, the two binary phenotypes |
| `15_fig_validation` | run, 1 s | Figure 1, framed as validation. Panel C shows the index recovers measured wing loading at R² = 0.93 |
| `16_rate_contrasts` | run, 6 s | Genome-wide dN, dS and omega contrasts on the index, with migration as a second orthogonal predictor |
| `17_power_sim` | run, 28 s | Detection power of the 42-pair design, computed analytically from the non-central t |
| `18_hyphaeon_prep` | run, 33 s | Phenotype table and gene manifest. 14,124 of 14,957 genes retained |
| `19_hyphaeon_run.py` | run, 26 h | The per-gene screen. Python, not R. See below |
| `20_hyphaeon_collect` | run | Parses the JSON, false-discovery control, base-composition check, non-neural comparator |
| `20b_contrast_screen` | run | **Primary per-gene analysis.** Sister-pair rate contrasts against the flight-index contrast. No phylogenetic null, so nothing to misspecify |
| `21_gsea` | run | GSEA on the continuous signed statistic, no threshold cut |
| `22_age_latency` | run, 2 s | Pair age and trait interaction, minimal-change baseline, latency simulation |
| `23_mrpmm` | run, 9.2 h per fit | Multi-response phylogenetic mixed model. Four specifications, see below |
| `24_fig_genomewide` | run | Figure 6, the genome-wide excess against a permutation null |
| `25_mrpmm_compare` | run | Compares the four MR-PMM fits and labels each correlation |
| `26_pelican_prep` | run, 3 min | Inputs for the per-site PELICAN screen. See `hpc/` |

### Not carried forward

Scripts 14 to 16 of the previous repository. The frequentist and Bayesian path analyses are
replaced by `23_mrpmm`, and their outputs are not included.

## The per-gene screen

`19_hyphaeon_run.py` **targets this Intel Mac specifically and says so at the top.**
PyTorch's last macOS x86_64 build is 2.2.2, which predates NumPy 2, so the whole stack must
be pinned to the NumPy 1 era before hyphaeon is installed. The recipe, the reasoning and a
startup environment check are in the script and in `HYPHAEON.md`. On Apple Silicon or Linux
none of this is needed.

It drives the library directly rather than calling the CLI per gene: analysis is 1.5 s but
process startup is 5.3 s, so 45,000 CLI invocations would waste roughly 65 hours.

Measured cost: **5.07 s per gene**, so about 20 hours per phenotype, and about 20 hours wall
for all three run as parallel processes. Permulations are free (0, 100 and 1000 take the
same time), so the default is 1000. `--resume` makes the job restartable.

HyphAeon's gene-level inference is RERconverge-style, but **it is not used for inference
here**. Its permulation null assumes pure Brownian motion and is exactly calibrated only
when the trait is Brownian; our flight index sits at lambda 0.85, which makes the test
conservative (whole-screen median p 0.786, zero genes at q < 0.05). Its gene summaries are
also not null-centred, so they cannot rank genes for GSEA. `HYPHAEON.md` documents the
calibration experiment. The primary per-gene analysis is `20b_contrast_screen`, which uses
sister-pair contrasts and therefore needs no phylogenetic null at all.

The RERconverge package itself is not used. It would supply no new null, only a different
statistic, and that statistic is recovered for free in `20_hyphaeon_collect` by
double-centring the log per-locus per-taxon rates in `dnds_branches.Rdata`, which is the
residual RERconverge models. 13,625 genes in minutes, from our own PAML rates.

## The multi-response model

`23_mrpmm` fits seven per-taxon responses with a phylogenetic random effect correlated
across responses. **Residual correlation is off.** It was retained at first, did not
converge, and is kept under the `_rescor` suffix as evidence.

Restoring it was reconsidered once outlier trimming had fixed the skew, and ruled out on
measurement. Trimming raises Pagel's lambda for ln dS from 1.025 to 1.079, for the flight
index from 0.960 to 1.026 and for GC from 0.912 to 1.004. The outliers were the main
departure from the tree, so removing them leaves data that look more Brownian and pushes
all three further onto the boundary where residual variance wants to be zero.

### The four specifications

| fit | flight response | taxa | reported correlations |
|---|---|---|---|
| `_norescor` | loading index | 83 | Rhat 1.047, ESS 79 |
| `_norescor_rz4` | loading index | 80 | **Rhat 1.011, ESS 543** |
| `_norescor_mig` | binary migratory | 83 | Rhat 1.012, ESS 335 |
| `_norescor_mig_rz4` | binary migratory | 80 | Rhat 1.018, ESS 209 |

Every parameter above Rhat 1.05 in every fit is a residual standard deviation or a term
sitting on ln dS, and none is a reported quantity. Four correlations are credible in all
four fits and agree to within 0.05: genome length with intron length at +0.68, dS with
omega at -0.59 (mechanical, since omega carries dS in its denominator), GC with omega at
-0.47 and GC with gene count at +0.37.

**Neither flight axis correlates with anything**, under either coding, trimmed or not.
Migration's phylogenetic standard deviation spans [0.08, 4.71], because its lambda is 0.00
and a phylogenetic covariance model has nothing to estimate from a trait with no
phylogenetic signal. Its variance lives at the tips, which is what the sister pairs read.

The first specification used ln dN and ln dS and would not converge (Rhat 1.99, 40
divergences, an hour without finishing). Two measured causes: ln dN and ln dS correlate at
**r = 0.970**, and Pagel's lambda was at or above 1 for three of four responses, meaning
residual variance wanted to be zero. Using **dS and omega** instead drops collinearity to
0.505 and gives omega lambda 0.00, so the residual term has real variance to estimate.

Responses, with measured phylogenetic signal:

| response | lambda |
|---|---|
| ln dS | 1.03 |
| ln omega | 0.00 |
| ln(mass / wing length²) | 0.96 |
| ln genome length | 0.60 |
| ln gene count | 0.31 |
| GC content | 0.91 |
| ln mean intron length | 0.00 |

Maximum pairwise |r| is 0.574. Dropped on evidence: ln dN (r = 0.97 with ln dS) and
single-exon gene count (r = 0.80 with gene count). BUSCO completeness is the only fixed
effect, since assembly quality does not evolve on the tree. Body mass is excluded because
it sits inside the flight index, which is a response.

The script reports skew, kurtosis, lambda and collinearity **before** fitting, and warns
rather than passing silently when Rhat exceeds 1.05. One response is flagged: the flight
index has skew 1.67 and kurtosis 7.08, driven by the wing-propelled divers and the ratites.
That is real biology, but a Gaussian response is an approximation there.

*Pinguinus* is excluded from this model only. Its assembly is 19 Mb of reads mapped to the
razorbill, so its genome length and gene count are mapping artefacts.

## Data

`data/Table_S1.xlsx` is the corrected trait table, with three sheets: the pair-level
Table S1 with AviList taxonomy and natural-log columns; **Table S2. Morphometrics**, 84
rows by 44 columns keyed by `Pair ID (analysis)` so it joins directly to `genomecomps`; and
**Table S3. Wing areas**, the 17-taxon validation set.

This replaces the previous repository's copy, which carried a body mass of 3560 g for
*Laterallus rogersi*, about ninety times too high, propagated from the *Pinguinus* cell.
*Pinguinus* body mass is 3561 g following Montgomerie and Birkhead (2025), replacing the
earlier PGLS estimate of 4750 g from Smith (2016).

Expensive intermediates are included rather than regenerated: `dnds_sisdifs_vol.Rdata`,
`dnds_branches.Rdata`, `p_values_per_locus.Rdata` and the Aitchison, genome-size and
permutation outputs.

`functions/` includes `trimCols.R` and `trimCodons.R`, which `02_align_basic_iq` sources
but which were absent from the previous repository.

**Codon alignments and gene trees are not in the repository.** They live at
`~/Desktop/Research/auk/postshrink_alignments` (29,914 files, about 14,957 genes, each
`*.fasta` with a matching `*.fasta.treefile`). Set `AUK_ALIGNMENTS` to point elsewhere.

### Matching our taxa to `Feng_Table_1.csv`

Two rules are needed, both established by inspection and both implemented in `23_mrpmm`:

1. Feng carries **four trinomials**, including `Aquila chrysaetos canadensis` and
   `Strix occidentalis caurina`. Matching on the full name silently loses those taxa.
   Compare the **first two words only**.
2. Feng's `Nannopterum harrisi` is `Phalacrocorax_harrisi` on our tree.

Together these recover 83 of 84. Only *Pinguinus* is absent, correctly, since it is not a
B10K genome.

## Reproducing this

Scripts run from `scripts/`, in numerical order, with the working directory set there.
Paths are relative: `../data/`, `../functions/`, `../figures/`, `../output/`.

Scripts 01 to 08 rebuild the alignments and rate estimates from the B10K data and take
days. They are included for completeness and for the record of how the inputs were made.
**Their expensive outputs are committed**, so scripts 09 onward run from this repository
alone, with the single exception of the codon alignments.

| stage | needs | typical cost |
|---|---|---|
| 01 to 08, alignment and rates | B10K c123 data, external tools below | days |
| 09 to 13, per-gene tests and genome size | committed `.Rdata` only | minutes |
| 14 to 22, traits, contrasts, screens | committed `.Rdata` only | minutes, except 19 |
| 19, the HyphAeon screen | codon alignments, hyphaeon | 26 h |
| 23, the MR-PMM | committed `.Rdata` only | 9.2 h per specification |
| 24 to 26, figures and PELICAN inputs | outputs of the above | minutes |

### External software

Carried over from the first version of this analysis, and needed only for scripts 01 to 08.

| | used by |
|---|---|
| MACSE v2.03 | 02, codon-aware alignment |
| IQ-TREE v1.7 | 02, gene trees |
| PAML | 04, dN and dS trees |
| EraBLE | 04, branch length estimation |
| BLASTN | 06, locus identification |

### R packages

`ape`, `phangorn`, `phytools`, `readxl`, `brms`, `jsonlite`, `lavaan`, `MASS`, `ggplot2`,
`gridExtra`, `clusterProfiler` 4.18.4, `fgsea` 1.36.2. R 4.5.2.

Gene sets are built from `data/goterms.csv` directly, 3,025 GO terms of which 389 fall
within the 10 to 500 size bounds, so no OrgDb is required and the annotation matches the
previous submission exactly.

### Helper functions

`functions/` holds `trimCols.R`, `trimCodons.R`, `trimLongTipBrs.R`, `rtt.cov.R`,
`pathnode.R` and `concatenate.R`. The first two were missing from the previous repository
although `02_align_basic_iq` sources them.

## Data provenance

| file | origin |
|---|---|
| `Table_S1.xlsx` | curated flight traits, AVONET morphometrics, AviList taxonomy |
| `dnds_sisdifs_vol.Rdata` | script 05, sister-pair rate ratios |
| `dnds_branches.Rdata` | script 04, per-locus per-taxon terminal branch rates |
| `p_values_per_locus.Rdata` | script 09, 14,785 by 24 per-gene binomial p-values |
| `Feng_Table_1.csv` | Feng et al. 2020, B10K assembly metrics |
| `goterms.csv` | GO annotation of the chicken reference loci |
| `dated.auk.tre` | the dated 84-taxon species tree |
| `OUT-0003_..._greatAuk.fa` | the *Pinguinus impennis* sequence, reads mapped to razorbill |
| `sheard_martin_traits.csv` | Sheard et al. hand-wing index and migration |
| `locus_codes.csv`, `taxonomy_codes.csv` | locus and taxon identifier maps |

**The codon alignments are not in this repository.** They are 29,914 files across roughly
14,957 genes, each `*.fasta` with a matching `*.fasta.treefile`, and they live at
`~/Desktop/Research/auk/postshrink_alignments`. Set `AUK_ALIGNMENTS` to point elsewhere.

## Deliberately not used

`RERconverge`, superseded by the in-house comparator in `20_hyphaeon_collect`, which adds
no new null. `MCMCglmm`, since brms serves. `ggtree`. HyphAeon is installed and the full
screen is retained, but it is not used for inference, for the reasons in `HYPHAEON.md`.

## Per-site screen

`hpc/` is a complete kit for running PELICAN on GenomeDK. It cannot be installed on the
machine this analysis was developed on, so the kit targets the cluster through Apptainer.
See `hpc/README.md`.
