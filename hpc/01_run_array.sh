#!/bin/bash
#SBATCH --account=DeiC-KU-L04
#SBATCH --job-name=pelican
#SBATCH --array=1-100
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/pelican_%A_%a.out
#
# The PELICAN screen, as a Slurm job array over gene chunks.
#
# Every gene is independent, so this scales almost linearly and coordinates nothing. Each
# array task takes a contiguous slice of the eligible gene list and writes into its own
# output directory, which 02_collect.R gathers afterwards.
#
# SET --account ABOVE AND --time FROM THE FIGURE 00_setup.sh PRINTED, THEN SUBMIT.
#
# Phenotype is chosen with PHENO, which defaults to flight:
#     sbatch --export=ALL,PHENO=flight 01_run_array.sh
#     sbatch --export=ALL,PHENO=loading 01_run_array.sh      # continuous, ~40 s/gene on 8 cores
# loading needs work/tree_loading.nhx from 01c_annotate_continuous.R first.
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"

PHENO=${PHENO:-flight}
source pelican_mode.txt                 # MODE, SECONDS_PER_GENE, PELICAN_CMD
PEL="$PELICAN_CMD"
# One alignment basename per line. Default is the flight-eligible set (>=70 taxa, >=12 of
# 15 reduced-flight tips); loading and migration use genes_min20.txt, the >=20-taxa rule of
# script 18. TAG names the output chunks so extra runs can sit beside earlier ones.
GENES=${GENES:-inputs/eligible_genes.txt}
NTASK=${SLURM_ARRAY_TASK_COUNT:-100}
TASK=${SLURM_ARRAY_TASK_ID:-1}
OUT=results/${PHENO}/${TAG:-chunk}_$(printf '%03d' "$TASK")
mkdir -p "$OUT" logs

# No pre-filter by default. The filters (multinomial for discrete, ANOVA for continuous)
# save time, but only sites passing p < 0.1 get the full test, so the tested p-values are
# selected and Benjamini-Hochberg over them is not valid. Set PREFILTER=0.1 to use them.
PREFILTER=${PREFILTER:-}
if [[ $PHENO == loading ]]; then
  SUB="continuous"; FILTER="--branch-trait average${PREFILTER:+ --anova-filter $PREFILTER}"
else
  SUB="discrete";   FILTER="${PREFILTER:+--multinomial-filter $PREFILTER}"
fi

# Contiguous slice of the gene list for this task.
TOTAL=$(wc -l < "$GENES")
PER=$(( (TOTAL + NTASK - 1) / NTASK ))
FROM=$(( (TASK - 1) * PER + 1 ))
TO=$(( TASK * PER )); (( TO > TOTAL )) && TO=$TOTAL
(( FROM > TOTAL )) && { echo "nothing to do for task $TASK"; exit 0; }
echo "task $TASK/$NTASK  genes $FROM..$TO of $TOTAL  mode $MODE  phenotype $PHENO"

SCRATCH=${TMPDIR:-/tmp}/pelican_$$_${TASK}
mkdir -p "$SCRATCH"
trap 'rm -rf "$SCRATCH"' EXIT

if [[ $MODE == BATCH ]]; then
  # One shared annotated tree, one invocation over a directory of alignments. PELICAN
  # requires that directory to contain ONLY alignments, hence the scratch copy. The copy is
  # uppercased because PELICAN silently reads lowercase bases as missing data, which turns
  # every site into "No observation" and still exits 0.
  mkdir -p "$SCRATCH/aln"
  sed -n "${FROM},${TO}p" "$GENES" | while read -r g; do
    awk '/^>/{print; next}{print toupper($0)}' "alignments/$g" > "$SCRATCH/aln/$g"
  done
  $PEL scan "$SUB" \
    --tree "work/tree_${PHENO}.nhx" --alignment "$SCRATCH/aln" \
    --alphabet nuc --exchangeability WAG $FILTER \
    --threads "${SLURM_CPUS_PER_TASK:-8}" -o "$OUT"
  echo "exit $?"
else
  # PERGENE: each gene carries its own taxon set, so each needs its own pruned, annotated
  # tree. PELICAN is a compiled binary, so the extra process per gene costs milliseconds.
  n=0; fail=0
  while read -r g; do
    base=${g%%.*}
    t=inputs/pruned_trees/${base}.nwk
    [[ -f $t ]] || { echo "SKIP $base, no pruned tree"; continue; }
    $PEL tree-annotation parsimony --tree "$t" \
         --annotation "inputs/pheno_${PHENO}.tsv" --clip \
         -o "$SCRATCH/${base}.nhx" >/dev/null 2>&1 \
      || { echo "SKIP $base, annotate failed"; fail=$((fail+1)); continue; }
    $PEL scan "$SUB" --tree "$SCRATCH/${base}.nhx" --alignment "alignments/$g" \
         --alphabet nuc --exchangeability WAG $FILTER \
         --threads "${SLURM_CPUS_PER_TASK:-8}" -o "$OUT/$base" >/dev/null 2>&1 \
      || { echo "SKIP $base, scan failed"; fail=$((fail+1)); continue; }
    n=$((n+1))
    (( n % 100 == 0 )) && echo "  ... $n genes done"
  done < <(sed -n "${FROM},${TO}p" "$GENES")
  echo "task $TASK finished: $n genes, $fail failures"
fi
