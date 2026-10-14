#!/bin/bash
#SBATCH --account=DeiC-KU-L04
#SBATCH --job-name=kegg_sets
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=../logs/%x_%j.out
#
# 21a_kegg_sets.sh : map each GALGAL_R locus to KEGG pathways for gallus gallus (gga).
#
# The loci carry B10K identifiers with no gene symbol or Entrez ID, so the chicken sequence
# in each alignment is translated and searched with DIAMOND against the RefSeq proteins of
# GRCg7b. The best hit gives an NCBI GeneID, which is the KEGG gga gene identifier.
# This replaces the KOBAS web server used in the previous submission, which annotated by
# the same kind of sequence search.
#
# Submit from hpc/kegg with `sbatch 21a_kegg_sets.sh`.
# Out: ../../data/galgal_refseq_map.csv, ../../data/kegg_sets.csv
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
export PATH=/faststorage/project/DeiC-KU-L04/shared_software/pixi/current/bin:$PATH
# The shared cache under shared_software/pixi holds files other members cannot read.
export PIXI_CACHE_DIR=/faststorage/project/DeiC-KU-L04/dduchene/.pixi_cache
# Bioconductor data packages (GenomeInfoDbData, GO.db) fetch their data in post-link scripts,
# which pixi skips by default, leaving clusterProfiler unloadable.
pixi config set --local run-post-link-scripts insecure
ALN=/faststorage/project/DeiC-KU-L04/dduchene/alignments_all
W=../work/kegg; mkdir -p "$W"
T=${SLURM_CPUS_PER_TASK:-2}
ASM=GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b
URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/$ASM

pixi install --manifest-path pixi.toml
MAN=$PWD/pixi.toml
run () { pixi run --manifest-path "$MAN" "$@"; }

echo "== chicken sequence from each alignment"
find "$ALN" -name 'GALGAL_R*.fasta' -print0 | xargs -0 awk '
  FNR == 1 { if (s != "") print ">" id "\n" s; s = ""; id = FILENAME; sub(/^.*\//, "", id); sub(/[.]c123.*$/, "", id); keep = 0 }
  /^>/     { keep = ($0 == ">GALGAL"); next }
  keep     { gsub(/[-!?.]/, ""); s = s $0 }
  END      { if (s != "") print ">" id "\n" s }' > "$W/galgal_nt.fa"
run seqkit translate -w 0 "$W/galgal_nt.fa" > "$W/galgal_aa.fa"
echo "loci with a chicken sequence: $(grep -c '>' "$W/galgal_aa.fa")"

echo "== RefSeq proteins and feature table"
[[ -s $W/gga.faa.gz ]] || curl -fsSL -o "$W/gga.faa.gz" "$URL/${ASM}_protein.faa.gz"
[[ -s $W/gga_features.txt.gz ]] || curl -fsSL -o "$W/gga_features.txt.gz" "$URL/${ASM}_feature_table.txt.gz"

echo "== DIAMOND"
run diamond makedb --in "$W/gga.faa.gz" -d "$W/gga" --threads "$T" --quiet
run diamond blastp -q "$W/galgal_aa.fa" -d "$W/gga" -o "$W/hits.tsv" --threads "$T" \
  --max-target-seqs 1 --evalue 1e-10 --sensitive --quiet \
  --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore

echo "== KEGG gga pathways"
curl -fsSL -o "$W/kegg_link.tsv" https://rest.kegg.jp/link/pathway/gga
curl -fsSL -o "$W/kegg_list.tsv" https://rest.kegg.jp/list/pathway/gga

cd ../../scripts
run Rscript 21a_kegg_sets.Rscript
