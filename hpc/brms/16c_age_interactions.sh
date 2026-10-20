#!/bin/bash
#SBATCH --account=DeiC-KU-L04
#SBATCH --job-name=age_interactions
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=../logs/%x_%j.out
#
# 16c_age_interactions.sh : run scripts/16c_age_interactions.Rscript (16 multi-response brms
# fits) in the shared R container. Submit from hpc/brms with `sbatch 16c_age_interactions.sh`.
# brms 2.22 is installed in R_LIBS_USER for the r-ver container.
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
T=/faststorage/project/DeiC-KU-L04/dduchene/tools
R="apptainer exec --env R_LIBS_USER=$T/Rlib /faststorage/project/DeiC-KU-L04/shared_software/containers/r-ver_4.4.1.sif Rscript"
cd ../../scripts
st=$(date +%s)
$R 16c_age_interactions.Rscript
echo "rc $? in $(( $(date +%s) - st )) s"
