#!/bin/bash
#SBATCH --account=DeiC-KU-L04
#SBATCH --job-name=traits_impute
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=../logs/%x_%j.out
#
# 14_traits.sh : run scripts/13b_generation_length.Rscript and scripts/14_traits_assemble.Rscript
# in the pixi environment of this folder (phytools needs GLPK, absent from the r-ver container).
# Submit from hpc/impute with `sbatch 14_traits.sh`.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
export PATH=/faststorage/project/DeiC-KU-L04/shared_software/pixi/current/bin:$PATH
export PIXI_CACHE_DIR=/faststorage/project/DeiC-KU-L04/dduchene/.pixi_cache
MAN=$PWD/pixi.toml
pixi install --manifest-path "$MAN"
cd ../../scripts
pixi run --manifest-path "$MAN" Rscript 13b_generation_length.Rscript
pixi run --manifest-path "$MAN" Rscript 14_traits_assemble.Rscript
