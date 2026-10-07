#!/bin/bash
#SBATCH --account=DeiC-KU-L04
#SBATCH --job-name=pelican_setup
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=logs/%x_%j.out
#
# 00_setup.sh : pull PELICAN, verify it runs, and decide which mode the screen uses.
#
# Submit from hpc/ with `sbatch 00_setup.sh`; nothing here runs on the login node. It is
# cheap and it resolves the one design question the screen depends on. In PERGENE mode it
# also writes the pruned trees, so the array can follow without another round trip.
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"

VER=1.0.8
SW=/faststorage/project/DeiC-KU-L04/shared_software/containers
IMG=$SW/pelican_${VER}.sif
RIMG=$SW/r-ver_4.4.1.sif
# Branch lengths must be substitutions per site. tree_codes.nwk is the dated tree in Myr,
# which inflates the test badly (see 00b_calibration.sh); tree_subst.nwk is it times 0.001.
TREE=inputs/tree_subst.nwk
PHENO=inputs/pheno_flight.tsv
THREADS=${SLURM_CPUS_PER_TASK:-2}
export APPTAINER_CACHEDIR=${TMPDIR:-/tmp}/apptainer_cache

say () { printf '\n=== %s ===\n' "$*"; }

# The alignments live outside the repository, in the user's folder; link rather than copy.
[[ -e alignments ]] || ln -s ../../alignments alignments
mkdir -p "$SW" work logs

say "1. Apptainer images"
[[ -f $IMG ]] || apptainer pull "$IMG" "docker://gitlab-registry.in2p3.fr/phoogle/pelican:${VER}" \
  || { echo "FATAL: could not pull the PELICAN image"; exit 1; }
[[ -f $RIMG ]] || apptainer pull "$RIMG" docker://rocker/r-ver:4.4.1 \
  || { echo "FATAL: could not pull the R image"; exit 1; }
PEL="apptainer exec $IMG pelican"
$PEL --version || { echo "FATAL: image will not run"; exit 1; }

# R with ape, for pruning and collection. rocker/r-ver installs CRAN binaries from Posit
# Package Manager, so this takes seconds and compiles nothing.
RLIB=$PWD/work/Rlib; mkdir -p "$RLIB"
RS="apptainer exec --env R_LIBS_USER=$RLIB $RIMG Rscript"
$RS -e 'if (!requireNamespace("ape", quietly=TRUE)) install.packages("ape", lib=Sys.getenv("R_LIBS_USER")); library(ape); cat("ape", format(packageVersion("ape")), "\n")' \
  || { echo "FATAL: could not install ape"; exit 1; }

say "PELICAN help, kept for reference"
for c in "tree-annotation parsimony" "scan discrete" "scan continuous"; do
  echo "--- pelican $c --help"; $PEL $c --help 2>&1
done

say "2. Annotate the full 84-tip tree with flight states"
# parsimony reconstructs ancestral discrete states at internal nodes from the leaf labels.
rm -f work/tree_flight.nhx    # PELICAN refuses to overwrite, and a stale tree is silent
$PEL tree-annotation parsimony \
  --tree "$TREE" --annotation "$PHENO" -o work/tree_flight.nhx \
  || { echo "FATAL: annotate-tree failed"; exit 1; }
echo "annotated tree written"

say "3. Smoke test A, a gene carrying ALL 84 taxa"
# If this fails, something is wrong with the inputs rather than with taxon coverage.
FULL=$(awk -F, 'NR>1 && $3==84 {gsub(/"/,"",$1); print $1; exit}' inputs/gene_manifest_eligible.csv)
FULL=alignments/$(basename "$FULL")
echo "gene: $(basename "$FULL")"
rm -rf work/smokeA
START=$(date +%s)
if $PEL scan discrete --tree work/tree_flight.nhx --alignment "$FULL" \
      --alphabet nuc --exchangeability WAG --threads "$THREADS" -o work/smokeA 2>work/smokeA.err; then
  echo "PASS"
else
  echo "FAIL. Inputs are wrong, not just taxon coverage. Error follows:"; cat work/smokeA.err; exit 1
fi
ELAPSED=$(( $(date +%s) - START ))

say "4. Smoke test B, THE DECIDING TEST: a gene missing some of the 84"
# 7087 of 9890 eligible genes have a unique taxon set, so if PELICAN insists that tree and
# alignment match, batching is impossible and every gene needs its own pruned tree.
PART=$(awk -F, 'NR>1 && $3>=70 && $3<=78 {gsub(/"/,"",$1); print $1; exit}' inputs/gene_manifest_eligible.csv)
PART=alignments/$(basename "$PART")
echo "gene: $(basename "$PART"), taxa: $(grep -c '^>' "$PART") of 84"
rm -rf work/smokeB
if $PEL scan discrete --tree work/tree_flight.nhx --alignment "$PART" \
      --alphabet nuc --exchangeability WAG --threads "$THREADS" -o work/smokeB 2>work/smokeB.err; then
  MODE=BATCH
  echo "PASS. PELICAN tolerates taxa absent from the alignment."
else
  MODE=PERGENE
  echo "FAIL. PELICAN requires tree and alignment to agree. Error follows:"
  head -5 work/smokeB.err
fi

say "5. Verdict"
echo "MODE=$MODE"                     | tee  pelican_mode.txt
echo "SECONDS_PER_GENE=$ELAPSED"      | tee -a pelican_mode.txt
echo "PELICAN_CMD=\"$PEL\""           | tee -a pelican_mode.txt
echo "RSCRIPT_CMD=\"$RS\""            | tee -a pelican_mode.txt
echo
echo "The 84-taxon gene took ${ELAPSED}s on $THREADS threads, unfiltered."
echo "9890 genes at that rate: about $(( ELAPSED * THREADS * 9890 / 3600 )) core-hours."

if [[ $MODE == PERGENE ]]; then
  say "6. Pruned trees for PERGENE mode"
  $RS 01b_prune_trees.R || { echo "FATAL: pruning failed"; exit 1; }
fi
echo "Set --time in 01_run_array.sh from the figure above, then: sbatch 01_run_array.sh"
