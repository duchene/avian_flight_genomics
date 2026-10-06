# Running PELICAN on GenomeDK

PELICAN tests, per protein site, whether amino acid preferences differ between foreground
and background lineages. It is the one resolution our analyses have not examined: a few
sites under strong selection in a few genes move neither genome-wide nor gene-level omega,
so our null on omega says nothing about them.

**Run this on GenomeDK, not locally.** The local Intel Mac cannot install it. Homebrew
refuses to build opam because the Command Line Tools predate Xcode 26.3, and Docker
Desktop 20.10.6 starts its backend but never opens a socket. Both are fixable only with
`sudo` and a system update. GenomeDK needs neither, because Apptainer pulls the authors'
prebuilt image.

## Why the cluster wins on more than convenience

| | this Mac | GenomeDK |
|---|---|---|
| install | blocked, needs a system update | `apptainer pull`, about two minutes |
| cores | 8 physical | a Slurm array across many nodes |
| 9,890 genes | days, single machine | hours, embarrassingly parallel |

The work is embarrassingly parallel. Every gene is independent, so a job array over gene
chunks scales almost linearly and there is nothing to coordinate.

## The one unknown, which step 00 resolves

Alignments carry a median of 77 of our 84 taxa, and **7,087 of the 9,890 eligible genes
have a taxon set no other gene shares**. So grouping genes by taxon set is useless, and
everything depends on one question:

> Does PELICAN accept a tree carrying taxa that are absent from the alignment?

- **If yes**, run in BATCH mode. One invocation per chunk of genes, one shared annotated
  84-tip tree, `--alignment` pointed at a directory. Fastest.
- **If no**, run in PERGENE mode. Prune the tree to each gene's taxa, annotate, scan one
  gene at a time. PELICAN is a compiled OCaml binary, so process startup is milliseconds
  rather than the 5.3 seconds that made per-gene invocation ruinous for HyphAeon. The cost
  is real but tolerable.

`00_setup.sh` answers the question and writes `MODE` to `pelican_mode.txt`. Do not skip it.

## Steps

```bash
# on the GenomeDK frontend, in your project directory
bash 00_setup.sh                      # pull image, smoke test, decide MODE, time one gene
sbatch 01_run_array.sh                # the screen itself, reads MODE
Rscript 02_collect.R                  # gather, FDR-correct, report
```

`00_setup.sh` prints a measured seconds-per-gene. Multiply by 9,890, divide by the array
width, and set `--time` in `01_run_array.sh` accordingly before submitting.

## Getting it onto the cluster

Clone the repository, then ship the alignments separately. `hpc/inputs/` is tracked, so
the clone already carries the relabelled tree, the three phenotype tables, the eligible
gene list and its manifest. Only the codon alignments need transferring, because they have
never been in the repository and should not be.

```bash
# here
bash hpc/make_transfer.sh                 # 113 Mb tarball of the 9,890 eligible alignments

# on the GenomeDK frontend
git clone <repo url> avian_flight_genomics
cd avian_flight_genomics/hpc
rsync -avP <this mac>:/tmp/auk_alignments.tar.gz .
mkdir -p alignments && tar -xzf auk_alignments.tar.gz -C alignments
bash 00_setup.sh
```

Cloning rather than shipping a self-contained payload means the cluster copy stays under
version control. Fixes made there come back by `git pull` instead of by hand, which
matters because the mode question in the next section may well force a change to
`01_run_array.sh` after the smoke test.

## Phenotypes, in the order worth running

1. **flight**, discrete, 15 tips across 14 pairs. Reduced flight ability, which is flight
   loss plus the curated limited-flight taxa. Flight loss alone is 7 taxa, too few to
   estimate a foreground profile.
2. **loading**, continuous. PELICAN 1.0.8 has `scan continuous`, so the index needs no
   binarisation and the categorical predictor the referees attacked stays buried. Costlier
   than the discrete test, and the authors recommend the `--anova-filter` pre-filter.
3. **migration**, discrete, 29 tips. Best powered and worst motivated: our pair-level
   result moves dN and dS together and leaves omega flat, which is a mutation-rate
   signature rather than the selection these methods detect.

## Caution on the great auk

*Pinguinus impennis* appears in only 40% of genes, and its assembly is reads mapped to the
razorbill, which is its own sister in this design. At the per-site level that is reference
bias pulling the pair together. Treat any hit resting on PINIMP-ALCTOR as suspect, and
check whether it survives dropping that pair.

## References

Tamuri, dos Reis, Hay and Goldstein (2009) PLoS Comput Biol 5(11), the original model.
Duchemin, Lanore, Veber and Boussau (2022) Mol Biol Evol, the benchmark and PELICAN.
