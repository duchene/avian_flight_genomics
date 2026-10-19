#!/bin/bash
#SBATCH --account=DeiC-KU-L04
#SBATCH --job-name=pair_gsea
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=../logs/%x_%j.out
#
# 27_pair_gsea.sh : run scripts/27_pair_gsea.Rscript, per-pair GSEA on omega and the
# tests of recurrence across pairs. Submit from hpc/kegg with `sbatch 27_pair_gsea.sh`.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
export PATH=/faststorage/project/DeiC-KU-L04/shared_software/pixi/current/bin:$PATH
# The shared cache under shared_software/pixi holds files other members cannot read.
export PIXI_CACHE_DIR=/faststorage/project/DeiC-KU-L04/dduchene/.pixi_cache
MAN=$PWD/pixi.toml
cd ../../scripts
pixi run --manifest-path "$MAN" Rscript 27_pair_gsea.Rscript
