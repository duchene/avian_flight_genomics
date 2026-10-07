#!/bin/bash
#SBATCH --account=DeiC-KU-L04
#SBATCH --job-name=pelican_checkcont
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source pelican_mode.txt
D=work/calib
rm -f work/tree_loading.nhx work/tree_loadperm.nhx
$RSCRIPT_CMD 01c_annotate_continuous.R loading || exit 1
paste <(tail -n +2 inputs/pheno_loading.tsv | cut -f1) <(tail -n +2 inputs/pheno_loading.tsv | cut -f2 | shuf --random-source=<(yes 1)) | sed '1i species\ttrait' > inputs/pheno_loadperm.tsv
$RSCRIPT_CMD 01c_annotate_continuous.R loadperm || exit 1
for lab in loading loadperm; do
  rm -rf $D/cont_$lab; s=$(date +%s)
  $PELICAN_CMD scan continuous --tree work/tree_$lab.nhx --alignment $D/aln --alphabet nuc \
    --exchangeability WAG --branch-trait average --threads 8 -o $D/cont_$lab >/dev/null 2>&1
  echo "$lab rc $? $(( $(date +%s)-s ))s"
done
rm -f inputs/pheno_loadperm.tsv
echo "labels   n  median  <0.05  <1e-3  <1e-6"
for lab in loading loadperm; do
  awk -F'\t' -v l=$lab 'NR==1{for(i=1;i<=NF;i++)if($i=="aagtr_pval")c=i} NR>1 && $c!="NA"{n++; v[n]=$c; a+=($c<0.05); b+=($c<1e-3); e+=($c<1e-6)} END{asort(v); printf "%-8s %d %.3f %.4f %.4f %.4f\n", l,n,v[int(n/2)+1],a/n,b/n,e/n}' $D/cont_$lab/all_sites.tsv
done
head -2 $D/cont_loading/all_sites.tsv; cut -f9 $D/cont_loading/all_sites.tsv | sort | uniq -c
