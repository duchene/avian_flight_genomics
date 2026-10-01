# avian_flight_genomics

Scripts and data for the analysis of avian flight habits and genome evolution, prepared
for submission to Molecular Biology and Evolution.

This repository replaces `avi_family_flight_genomics`. It carries forward the alignment,
dating and sister-pair pipeline unchanged, and replaces the categorical flight predictor,
the gene-set enrichment approach and the Bayesian path analysis. Nothing superseded has
been copied across, so a release from here contains no mixture of old and new material.

## Design

42 sister pairs, 84 taxa, roughly 14,800 protein-coding genes. All analyses are restricted
to these 84 taxa, which are the species with curated flight traits. No analysis uses the
wider 363-taxon B10K set.

The primary flight predictor is **ln(mass / wing length²)**, a wing-loading index,
replacing the four-level flight-habit categories used previously. Wing loading proper is
mass over wing area, and published wing areas exist for only 8 of the 84 taxa. Wing length
is available for all 84, from AVONET for 83 and from Livezey (1988, Table 1) for
*Pinguinus impennis*. Wing area scales as wing length to the power 1.9 to 2.4 where both
are known, so the index has the dimensions of wing loading and tracks measured loading
with R² = 0.92 across 17 g to 4664 g and six flight styles.

Two further axes enter as binary predictors: gain of migration, and flight loss. The
wing-loading index is statistically orthogonal to migration, so both can be fitted without
collinearity.

All logarithms are natural, matching `masdif` in `genome_flight_age_mass_data.Rdata`.

## Pipeline

Scripts run from `scripts/`, with paths `../data/`, `../functions/`, `../figures/` and
`../output/`.

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
| `10_aitchison_variance` | Base composition heterogeneity per gene |
| `11_base_composition_vs_rates` | Aitchison variance against binomial p-values |
| `12_ks_permutation_envelopes` | Permutation null envelopes for the KS tests |
| `13_genome_size` | Genome metrics and assembly-quality correction |

### New for this submission

| | |
|---|---|
| `14_traits_assemble` | **Written.** Builds the per-taxon trait table, the wing-loading index and its within-pair contrasts, the two binary phenotypes, and mean GC per taxon. GC variance across genes is computed only if `data/postshrink_alignments/` is present |
| `15_fig_validation` | **Written.** Figure 1. Index against primary flight style and against flight ability, the 8 taxa with measured wing area (R² = 0.93), and a ladderised tree coloured by flight style. Framed as validation of the index |
| `16_rate_contrasts` | Genome-wide dN, dS and ω contrasts against the index, with migration as a second orthogonal predictor |
| `17_power_sim` | Detection power of the 42-pair design, run before the per-gene screen |
| `18_hyphaeon_prep` | Write alignments, trees and phenotype files for three analyses: continuous index, binary flight loss, binary migration |
| `19_hyphaeon_run.sh` | Run `hyphaeon phenotype` across genes for each phenotype |
| `20_hyphaeon_collect` | Collect per-gene statistics, false-discovery control, and run RERconverge on the same phenotypes as the peer-reviewed parallel |
| `21_gsea` | clusterProfiler GSEA on the continuous per-gene statistic, no threshold cut |
| `22_age_latency` | Pair age and trait interaction, minimal-change baseline, and transition-latency simulation |
| `23_mrpmm` | Multi-response phylogenetic mixed model across the 84-taxon phylogeny |
| `24_figures` | Publication figures, one named function per figure |

### Not carried forward

Scripts 14 to 16 of the previous repository. The frequentist and Bayesian path analyses
are replaced by the multi-response phylogenetic mixed model in `23_mrpmm`, and their
outputs (`bayesian_path_results.Rdata`, `bayesian_path_sensitivity.Rdata`,
`bayesian_model_comparison*.csv`) are not included.

## Data

`data/Table_S1.xlsx` is the corrected trait table. It carries three sheets.

- **Table S1**, pair level, AviList taxonomy, natural-log columns
- **Table S2. Morphometrics**, 84 rows by 44 columns, one row per taxon, keyed by
  `Pair ID (analysis)` so that it joins directly to `genomecomps`
- **Table S3. Wing areas**, the 17-taxon wing-area validation set

This replaces the copy in the previous repository, which carried a body mass of 3560 g for
*Laterallus rogersi*, about ninety times too high, propagated from the *Pinguinus* cell.
*Pinguinus* body mass is 3561 g following Montgomerie and Birkhead (2025), replacing the
earlier PGLS estimate of 4750 g from Smith (2016).

Expensive intermediates are included rather than regenerated: `dnds_sisdifs_vol.Rdata`,
`dnds_branches.Rdata`, `p_values_per_locus.Rdata` and the Aitchison, genome-size and
permutation outputs.

`functions/` now includes `trimCols.R` and `trimCodons.R`, which `02_align_basic_iq`
sources but which were absent from the previous repository.

## Status

`STATUS.md` is the working handoff note: decisions already taken, numbers already
established, traps, and the next action. Read it before resuming work.

Scripts 01 to 13 are carried over and run. Scripts 14 and 15 are written and tested.
Scripts 16 to 24 are specified in `../PLAN_next_submission.md` but not yet written.

Script 14 writes `output/traits_per_taxon.csv` (84 rows) and `output/traits_contrasts.csv`
(42 rows). Every downstream script reads those two files rather than the spreadsheet.

Not installed here: `clusterProfiler` (needed by 21), `RERconverge` (needed by 20),
`ggtree`, `MCMCglmm`. Script 23 uses `brms`, which is installed.

## Open questions before the remaining scripts are written

1. Does `hyphaeon phenotype` require a phenotype value for every tip in an alignment, or
   only for the focal taxa? This determines how `18_hyphaeon_prep` prunes.
2. Is the HyphAeon per-gene output a rate-association statistic, as in RERconverge, or a
   diversifying-selection statistic, as in MEME? The preprint does not say, and the
   reading of a genome-wide excess as relaxed selection depends on the answer.
3. Runtime and batching for the full per-gene screen across three phenotypes.
4. Whether each gene analysis uses its own gene tree or the dated species tree.

Ancient-DNA damage control for the great auk has to be completed before `19_hyphaeon_run`,
because cytosine deamination produces C-to-T and G-to-A miscalls that bias dN and dS
directly, and *Pinguinus* enters the codon-level analyses.
