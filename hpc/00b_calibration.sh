#!/bin/bash
#SBATCH --account=DeiC-KU-L04
#SBATCH --job-name=pelican_calib
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
#
# 00b_calibration.sh : is the PELICAN test calibrated on our tree?
#
# tree_codes.nwk is the dated tree, so its branch lengths are millions of years, while
# PELICAN's model reads them as expected substitutions per site. The first smoke test gave
# p of 1e-16 at sites whose amino acid mix is the same in both groups. This job scans 20
# random genes with the tree rescaled by several factors, each under the real flight labels
# and under labels permuted across tips. A calibrated test gives roughly uniform p-values
# under permutation, so the permuted column is the one that decides the scale.
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source pelican_mode.txt
PEL="$PELICAN_CMD"; RS="$RSCRIPT_CMD"
SCALES=${SCALES:-"1 0.01 0.001 0.0003"}
D=work/calib; rm -rf "$D"; mkdir -p "$D/aln"

# PELICAN silently treats lowercase bases as missing, so every alignment is uppercased.
shuf -n 20 --random-source=<(yes) inputs/eligible_genes.txt | while read -r g; do
  awk '/^>/{print; next}{print toupper($0)}' "alignments/$g" > "$D/aln/$g"
done

# Permuted labels: same tips, same counts, shuffled states.
paste <(tail -n +2 inputs/pheno_flight.tsv | cut -f1) \
      <(tail -n +2 inputs/pheno_flight.tsv | cut -f2 | shuf --random-source=<(yes 1)) \
  | sed '1i species\ttrait' > "$D/pheno_perm.tsv"

for s in $SCALES; do
  $RS -e "library(ape); t <- read.tree('inputs/tree_codes.nwk'); t\$edge.length <- t\$edge.length * $s; write.tree(t, '$D/tree_s$s.nwk')"
  for lab in real perm; do
    ph=inputs/pheno_flight.tsv; [[ $lab == perm ]] && ph=$D/pheno_perm.tsv
    $PEL tree-annotation parsimony --tree "$D/tree_s$s.nwk" --annotation "$ph" -o "$D/tree_s${s}_$lab.nhx"
    START=$(date +%s)
    $PEL scan discrete --tree "$D/tree_s${s}_$lab.nhx" --alignment "$D/aln" \
      --alphabet nuc --exchangeability WAG --threads "${SLURM_CPUS_PER_TASK:-8}" \
      -o "$D/out_s${s}_$lab" >/dev/null 2>&1
    echo "scale $s $lab: rc $? in $(( $(date +%s) - START ))s"
  done
done

echo; echo "scale labels  n_tested  median_p  frac<0.05  frac<1e-3  frac<1e-6"
for s in $SCALES; do for lab in real perm; do
  awk -F'\t' -v s="$s" -v l="$lab" 'NR>1 && $6!="NA"{n++; v[n]=$6; a+=($6<0.05); b+=($6<1e-3); c+=($6<1e-6)}
    END{asort(v); printf "%-6s %-5s %9d %9.3f %10.3f %10.4f %10.4f\n", s, l, n, v[int(n/2)+1], a/n, b/n, c/n}' \
    "$D/out_s${s}_$lab/all_sites.tsv"
done; done
