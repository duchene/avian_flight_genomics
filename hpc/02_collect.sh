#!/bin/bash
#SBATCH --account=DeiC-KU-L04
#SBATCH --job-name=pelican_collect
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%j.out
#
# Run 02_collect.R inside the R container. Usage: sbatch 02_collect.sh flight
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source pelican_mode.txt
$RSCRIPT_CMD 02_collect.R "${1:-flight}"
