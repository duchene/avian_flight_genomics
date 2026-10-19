# PELICAN on GenomeDK

PELICAN tests, per protein site, whether amino acid preferences differ between conditions
or covary with a continuous trait. It is the one resolution the rest of the analysis does
not examine: a few sites under strong selection in a few genes move neither genome-wide nor
gene-level omega.

**Result: null for all three phenotypes.** No site and no gene reaches q < 0.10. Per-gene
tables are in `results/summary/`.

| | reduced flight (15 tips) | loading (continuous) | migration (29 tips) |
|---|---|---|---|
| genes | 9,890 eligible | 14,124 (>= 20 taxa) | 14,124 (>= 20 taxa) |
| variable sites tested | 1.90M | 2.79M | 2.79M |
| sites at BH q < 0.10 | 0 | 0 | 0 |
| genes at q < 0.10 (Sidak within, BH across) | 0 | 0, best q 0.129 | 0 |
| smallest site p | 1.2e-5 | 2.5e-7 | 8.8e-6 |
| median p, fraction < 0.05 | 0.57, 0.034 | 0.45, 0.036 | 0.35, 0.024 |

The near miss for loading, GALGAL_R08348 site 64, carries Asn in most low-loading taxa and
Ser in nearly all high-loading ones, but the Asn taxa are largely passerines and their
allies, so it reads as a single clade-level change.

Every step is a Slurm job on the `DeiC-KU-L04` account. Nothing runs on the login node.

## Steps

```bash
cd hpc
sbatch 00_setup.sh            # pull images, install ape, annotate flight tree, smoke tests
sbatch 00b_calibration.sh     # optional: test calibration under permuted labels
# loading only: write work/tree_loading.nhx (Brownian ancestral values)
srun --account=DeiC-KU-L04 --mem=2G --time=0:10:00 bash -c \
  'source pelican_mode.txt; $RSCRIPT_CMD 01c_annotate_continuous.R loading'
sbatch --export=ALL,PHENO=flight 01_run_array.sh
sbatch --export=ALL,PHENO=migration,GENES=inputs/genes_min20.txt 01_run_array.sh
sbatch --export=ALL,PHENO=loading,GENES=inputs/genes_min20.txt   01_run_array.sh
sbatch 02_collect.sh flight   # likewise migration, loading
```

The migration tree (`work/tree_migration.nhx`) is annotated like the flight tree, with
`pelican tree-annotation parsimony` on `inputs/tree_subst.nwk`.

## Inputs

| file | what |
|---|---|
| `inputs/tree_codes.nwk` | dated 84-taxon tree, Myr, six-letter codes. Used only for ancestral loading |
| `inputs/tree_subst.nwk` | the same tree times 0.001, substitutions per site. **What PELICAN reads** |
| `inputs/pheno_{flight,migration,loading}.tsv` | header `species<TAB>trait` |
| `inputs/eligible_genes.txt` | 9,890 genes, >= 70 taxa and >= 12 of 15 reduced-flight tips. Flight |
| `inputs/genes_min20.txt` | 14,124 genes, >= 20 taxa, the script 18 rule. Loading and migration |
| `alignments` | symlink to `dduchene/alignments_all`, all 14,957 `*_NT.fasta` |

Images live in `/faststorage/project/DeiC-KU-L04/shared_software/containers/`: PELICAN 1.0.8
and rocker R 4.4.1, with `ape` in `work/Rlib`. `pelican_mode.txt`, written by setup, holds
the two container commands every script sources.

## Choices, and the evidence for each

- **No pre-filter.** PELICAN's multinomial and ANOVA filters pass only sites already at
  p < 0.1 on a label-dependent screen, so BH over the tested sites is invalid. On flight the
  filter manufactured 22,911 sites at q < 0.10 that vanish without it. `PREFILTER=0.1`
  restores it if speed ever matters more than valid q-values.
- **Invariant sites dropped before BH** (`02_collect.R`). Invariance does not depend on the
  labels, so this is valid independent filtering and only removes sites that cannot be
  significant.
- **Branch lengths in substitutions.** On the dated tree, permuted labels put 1.5% of sites
  at p < 1e-6. Rescaled by 0.001 the test behaves, identically from x0.01 to x0.0003.
- **Loading's internal nodes** are Brownian-motion ML ancestral values from the dated tree,
  and branches take the mean of their ends (`--branch-trait average`). Lambda for loading is
  0.85, but a lambda-scaled reconstruction correlates with these at r = 0.999 (`checks/`).
- **BATCH mode.** PELICAN accepts the 84-tip tree for a gene missing taxa, so one shared tree
  serves every gene. `01b_prune_trees.R` and the PERGENE branch are unused.

## Silent traps in PELICAN 1.0.8

Each of these produced plausible output or a zero exit code, not an error.

- Lowercase bases are read as missing. Every site becomes "No observation" and the run
  exits 0. `01_run_array.sh` uppercases each alignment into scratch.
- The annotation header must be `trait`, though the help text says `traits`.
- The subcommand is `tree-annotation`, not `annotate-tree`.
- It will not overwrite an existing output file, so a stale annotated tree survives.

## Cost

Discrete scan about 2.3 s per gene on 8 cores, continuous about 40 s. The full set of five
arrays, setup and checks came to about 720 core-hours.

## Storage

Raw per-chunk output is `results/<phenotype>/*/all_sites.tsv(.gz)`. PELICAN's `store.pag`
resume database is not needed once a run has finished.

## References

Tamuri, dos Reis, Hay and Goldstein (2009) PLoS Comput Biol 5(11), the original model.
Duchemin, Lanore, Veber and Boussau (2022) Mol Biol Evol, the benchmark and PELICAN.
