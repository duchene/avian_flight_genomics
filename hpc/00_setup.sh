#!/bin/bash
# 00_setup.sh : pull PELICAN, verify it runs, and decide which mode the screen uses.
#
# Run this on the GenomeDK frontend before submitting anything. It is cheap and it
# resolves the one design question the screen depends on.
set -uo pipefail

VER=1.0.8
IMG=pelican_${VER}.sif
TREE=inputs/tree_codes.nwk
PHENO=inputs/pheno_flight.tsv
ALN=alignments

say () { printf '\n=== %s ===\n' "$*"; }

say "1. Apptainer image"
if [[ -f $IMG ]]; then
  echo "already present: $IMG"
else
  # GenomeDK ships Singularity/Apptainer and can execute Docker images directly.
  apptainer pull "$IMG" "docker://gitlab-registry.in2p3.fr/phoogle/pelican:${VER}" \
    || singularity pull "$IMG" "docker://gitlab-registry.in2p3.fr/phoogle/pelican:${VER}" \
    || { echo "FATAL: could not pull the image"; exit 1; }
fi
PEL="apptainer exec $IMG pelican"
command -v apptainer >/dev/null || PEL="singularity exec $IMG pelican"
$PEL --version || { echo "FATAL: image will not run"; exit 1; }

say "2. Annotate the full 84-tip tree with flight states"
# parsimony reconstructs ancestral discrete states at internal nodes from the leaf labels.
mkdir -p work
$PEL annotate-tree parsimony \
  --tree "$TREE" --annotation "$PHENO" -o work/tree_flight.nhx \
  || { echo "FATAL: annotate-tree failed"; exit 1; }
echo "annotated tree written"

say "3. Smoke test A, a gene carrying ALL 84 taxa"
# If this fails, something is wrong with the inputs rather than with taxon coverage.
FULL=$(awk -F, 'NR>1 && $3==84 {print $1; exit}' inputs/gene_manifest_eligible.csv)
FULL=alignments/$(basename "$FULL")
echo "gene: $(basename "$FULL")"
rm -rf work/smokeA
if $PEL scan discrete --tree work/tree_flight.nhx --alignment "$FULL" \
      --alphabet nuc --exchangeability WAG --threads 2 -o work/smokeA 2>work/smokeA.err; then
  echo "PASS"
else
  echo "FAIL. Inputs are wrong, not just taxon coverage. Error follows:"; cat work/smokeA.err; exit 1
fi

say "4. Smoke test B, THE DECIDING TEST: a gene missing some of the 84"
# 7087 of 9890 eligible genes have a unique taxon set, so if PELICAN insists that tree and
# alignment match, batching is impossible and every gene needs its own pruned tree.
PART=$(awk -F, 'NR>1 && $3>=70 && $3<=78 {print $1; exit}' inputs/gene_manifest_eligible.csv)
PART=alignments/$(basename "$PART")
echo "gene: $(basename "$PART"), taxa: $(grep -c '^>' "$PART") of 84"
rm -rf work/smokeB
START=$(date +%s)
if $PEL scan discrete --tree work/tree_flight.nhx --alignment "$PART" \
      --alphabet nuc --exchangeability WAG --threads 2 -o work/smokeB 2>work/smokeB.err; then
  MODE=BATCH
  echo "PASS. PELICAN tolerates taxa absent from the alignment."
else
  MODE=PERGENE
  echo "FAIL. PELICAN requires tree and alignment to agree. Error follows:"
  head -5 work/smokeB.err
fi
ELAPSED=$(( $(date +%s) - START ))

say "5. Verdict"
echo "MODE=$MODE"      | tee  pelican_mode.txt
echo "SECONDS_PER_GENE=$ELAPSED" | tee -a pelican_mode.txt
echo "PELICAN_CMD=$PEL" | tee -a pelican_mode.txt
echo
echo "One gene took ${ELAPSED}s on 2 threads."
echo "9890 genes single-threaded would be about $(( ELAPSED * 9890 / 3600 )) core-hours."
echo "Across a 100-task array that is roughly $(( ELAPSED * 9890 / 3600 / 100 + 1 )) hours of wall time."
echo
if [[ $MODE == PERGENE ]]; then
  echo "PERGENE mode needs a pruned tree per gene. Generate them with:"
  echo "   Rscript 01b_prune_trees.R      # needs R with ape, e.g. conda install -c conda-forge r-ape"
fi
echo "Set --time in 01_run_array.sh from the figure above, then: sbatch 01_run_array.sh"
