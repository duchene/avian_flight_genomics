#!/bin/bash
#SBATCH --account=DeiC-KU-L04
#SBATCH --job-name=mrpmm
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=72:00:00
#SBATCH --output=../logs/%x_%j.out
#
# mrpmm.sh : run scripts/23_mrpmm.Rscript for one specification, in the shared R container.
# Usage, from hpc/brms:
#   sbatch --job-name=mrpmm_load_gen --export=ALL,MRPMM_FLIGHT=loading,MRPMM_TRIM=none,MRPMM_COVAR=gen mrpmm.sh
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
T=/faststorage/project/DeiC-KU-L04/dduchene/tools
R="apptainer exec --env R_LIBS_USER=$T/Rlib,MRPMM_FLIGHT=$MRPMM_FLIGHT,MRPMM_TRIM=$MRPMM_TRIM,MRPMM_COVAR=$MRPMM_COVAR /faststorage/project/DeiC-KU-L04/shared_software/containers/r-ver_4.4.1.sif Rscript"
cd ../../scripts
st=$(date +%s)
$R 23_mrpmm.Rscript
rc=$?; echo "rc $rc in $(( $(date +%s) - st )) s"; exit $rc
