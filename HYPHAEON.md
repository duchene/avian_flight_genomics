# HyphAeon: environment, interface and confirmed behaviour

Everything below was established by installing hyphaeon 0.2.0 and running it to completion
on a real gene from this project. It is not inferred from the preprint.

## The key finding

**HyphAeon `phenotype` is a RERconverge-style trait-rate association at the gene level.**
The `--permulations` flag is documented in the CLI as

> Number of Brownian motion phylogenetic permulations for empirical p-values
> (RERconverge null model; default: 0 / parametric)

So the gene-level p-value is an empirical permulation p-value under a Brownian-motion null,
not a diversifying-selection test. The reading of a genome-wide or widespread signal as
relaxed selection is therefore available, and no change of plan is needed.

Per site it returns both a transformer LRT (`hyphaeon_lrt`) and a trait association
(`association_rho`), so the site-level statistic is a hybrid of the two. The gene-level
number is the permulation one.

## Installation on an Intel Mac

This machine is x86_64. PyTorch's last macOS x86_64 build is **2.2.2**, which predates
NumPy 2. HyphAeon's default resolution installs NumPy 2 and the run then dies at
`tensor.numpy()` with "Numpy is not available". Downgrading NumPy afterwards segfaults,
because scipy's wheel is built against NumPy 2.

**The whole stack must be pinned to the NumPy 1 era before hyphaeon is installed:**

```sh
python3 -m venv hvenv
hvenv/bin/pip install --upgrade pip
hvenv/bin/pip install "numpy==1.26.4" "scipy==1.11.4" "pandas==2.1.4" \
                      "matplotlib==3.8.4" "torch==2.2.2"
hvenv/bin/pip install --no-deps hyphaeon hyphaeon-core
hvenv/bin/pip install biopython tn93 safetensors huggingface-hub tqdm pyyaml click truststore
```

Verify with:

```sh
hvenv/bin/python -c "import torch, aeon_core.model, hyphaeon; print(torch.zeros(2).numpy())"
```

On Apple Silicon or Linux, a plain `pip install hyphaeon` should work, since torch is
current there. The pinning above is specific to this machine.

Model weights (about 8 MB) download from Hugging Face on first use and cache in
`~/.cache/huggingface/`. No authentication needed.

**Debugging trap.** The `exceptiongroup` backport in this dependency tree **segfaults
inside its own exception hook** while formatting errors, so every real error message is
replaced by a traceback from `exceptiongroup/_formatting.py` and then a crash. To see a
real error, reset the hook after importing it:

```python
import exceptiongroup, sys
sys.excepthook = sys.__excepthook__
```

## Interface

```sh
hyphaeon phenotype \
  -a  <in-frame codon FASTA or NEXUS> \
  -t  <Newick or NEXUS tree>            # optional; --no-tree/--use-tn93 uses TN93 instead
  -pf <phenotype CSV or TSV> \
  -tc <trait column name> \
  -sc <species column name> \
  --continuous                          # omit for binary; or use -fg/-bg instead of -pf
  --permulations <N>                    # RERconverge null; 0 means parametric
  --alpha <FDR threshold>               # default 0.05
  --min-taxa <N>                        # default 4
  --cpu \
  -o <output.json>  -c <output.csv>
```

Our existing files work unchanged. Alignments and per-gene trees are at
`~/Desktop/Research/auk/postshrink_alignments` (29,914 files, so about 14,957 genes, each
`*.fasta` with a matching `*.fasta.treefile`). Taxon labels are the six-letter codes, which
match the `code` column of `output/traits_per_taxon.csv`. On the test gene it reported
**100% taxon concordance, 76 shared taxa**, with no reformatting.

Identical sequences are collapsed before analysis (76 taxa became 75 haplotypes on the
test gene), which is worth noting when reconciling taxon counts.

## Output schema

`-c` writes site-level rows only. **Gene-level statistics are in the `-o` JSON**, so the
collector must parse JSON, not CSV.

Gene level: `gene_p_value_perm`, `spectral_energy`, `norm_spectral_ratio`, `max_assoc`,
`p_evd_length_adjusted`, `score_track_a`, `score_track_b`, `dual_track_composite`,
`compact_pars_signature`, `permulations_count`, `significant_sites_count`,
`coselection_pairs_count`, `trait_sectors_count`, `taxa_count`, `codon_count`.

Site level, in `sites` and in the CSV: `site`, `ref_aa`, `derived_aa`, `hyphaeon_lrt`,
`p_lrt`, `attribution_norm`, `fg_mean_attn`, `bg_mean_attn`, `association_rho`, `p_value`,
`p_assoc`, `p_assoc_parametric`, `p_assoc_perm`, `score`, `foreground_freq_pct`,
`background_freq_pct`, `q_value`.

Note that `association_rho` is signed, while the gene-level `max_assoc` is a magnitude. A
signed gene-level statistic for GSEA has to be built from the sites.

## Timing and the runner

Measured on the test gene: 76 codons, 75 haplotypes, 100 permulations.

- analysis itself: **1.5 s** (about 105 codons/s for attribution)
- wall clock per process: **5.3 s**

Process startup and the torch import dominate, so **the runner must load the model once
and loop over genes in a single process** rather than invoking the CLI per gene. Calling
the CLI 45,000 times would cost roughly 65 hours of pure overhead.

At 1.5 s per gene, 14,957 genes across 3 phenotypes is about 19 hours on one core and a
few hours parallelised.

## Worked example

```sh
hvenv/bin/hyphaeon phenotype \
  -a  postshrink_alignments/GALGAL_R00003.c123.linsi_macse_NT.fasta \
  -t  postshrink_alignments/GALGAL_R00003.c123.linsi_macse_NT.fasta.treefile \
  -pf pheno.csv -tc loading_index -sc species \
  --continuous --permulations 100 --cpu \
  -o out.json -c out.csv
```

returned, in 1.5 s: taxa 75, codon sites 76, spectral energy 0.1075, normalised spectral
ratio 0.1280, gene empirical p = 0.713 from 100 permulations, 0 sites significant at
FDR 0.05.
