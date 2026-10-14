#!/bin/bash
#SBATCH --account=DeiC-KU-L04
#SBATCH --job-name=gsea_go_kegg
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=../logs/%x_%j.out
#
# 21b_gsea.sh : run scripts/21_gsea.Rscript with GO and with KEGG gene sets.
# Submit from hpc/kegg after 21a_kegg_sets.sh, with `sbatch 21b_gsea.sh`.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
export PATH=/faststorage/project/DeiC-KU-L04/shared_software/pixi/current/bin:$PATH
# The shared cache under shared_software/pixi holds files other members cannot read.
export PIXI_CACHE_DIR=/faststorage/project/DeiC-KU-L04/dduchene/.pixi_cache
# Bioconductor data packages (GenomeInfoDbData, GO.db) fetch their data in post-link scripts,
# which pixi skips by default, leaving clusterProfiler unloadable.
pixi config set --local run-post-link-scripts insecure
MAN=$PWD/pixi.toml
cd ../../scripts
for s in GO KEGG; do
  echo; echo "########## $s ##########"
  GSEA_SETS=$s pixi run --manifest-path "$MAN" Rscript 21_gsea.Rscript
done
