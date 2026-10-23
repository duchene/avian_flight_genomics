#!/bin/bash
#SBATCH --account=DeiC-KU-L04
#SBATCH --job-name=sens
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=../logs/%x_%j.out
#
# sens.sh : run one regression script in one sensitivity mode, in the shared R container.
# Usage, from hpc/brms:
#   sbatch --job-name=s16_gen --export=ALL,SCRIPT=16_rate_contrasts,SENS=gen sens.sh
# SENS=gen uses generation length in place of body mass; SENS=gc adds GC as a response.
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
T=/faststorage/project/DeiC-KU-L04/dduchene/tools
R="apptainer exec --env R_LIBS_USER=$T/Rlib,SENS=$SENS /faststorage/project/DeiC-KU-L04/shared_software/containers/r-ver_4.4.1.sif Rscript"
cd ../../scripts
st=$(date +%s)
$R "$SCRIPT.Rscript"
rc=$?; echo "rc $rc in $(( $(date +%s) - st )) s"; exit $rc
