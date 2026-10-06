#!/bin/bash
# make_transfer.sh : build the alignment payload for GenomeDK. Run locally, from the repo root.
#
# Everything else reaches the cluster by `git clone`, because hpc/inputs/ is tracked. Only
# the 9,890 codon alignments need shipping, since they have never been in the repository
# and should not be: they are 316 Mb raw across 14,957 files in their full form.
set -euo pipefail

OUT=${1:-/tmp/auk_alignments.tar.gz}
ALN=${AUK_ALIGNMENTS:-$HOME/Desktop/Research/auk/postshrink_alignments}
LIST=hpc/inputs/eligible_genes.txt
[[ -f $LIST ]] || { echo "FATAL: $LIST missing. Run scripts/26_pelican_prep.Rscript first."; exit 1; }

N=$(wc -l < "$LIST")
echo "packing $N eligible alignments from $ALN"
tar -czf "$OUT" -C "$ALN" -T "$LIST"
echo "payload: $OUT  $(du -h "$OUT" | cut -f1)"
cat <<MSG

On GenomeDK:

  git clone <your repo url> avian_flight_genomics
  cd avian_flight_genomics/hpc
  rsync -avP <this mac>:$OUT .          # or scp from here
  mkdir -p alignments && tar -xzf $(basename "$OUT") -C alignments
  bash 00_setup.sh

MSG
