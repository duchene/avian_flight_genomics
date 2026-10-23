#!/bin/bash
#SBATCH --account=DeiC-KU-L04
#SBATCH --job-name=gc3
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=../logs/%x_%j.out
#
# gc3.sh : GC3 per gene and taxon (scripts/19b_gc3.Rscript), its per-gene screen
# (scripts/20b_contrast_screen.Rscript with GC3=TRUE), and its genome-wide excess against
# the permutation null (auk_ms figure_code/fig3_genomewide.R with GC3=TRUE).
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
D=/faststorage/project/DeiC-KU-L04/dduchene
T=$D/tools
C=/faststorage/project/DeiC-KU-L04/shared_software/containers/r-ver_4.4.1.sif
ENV="R_LIBS_USER=$T/Rlib,AUK_ALIGNMENTS=$D/alignments_all,GC3=TRUE,AFG=$D/avian_flight_genomics"
cd ../../scripts
apptainer exec --env "$ENV" $C Rscript 19b_gc3.Rscript
apptainer exec --env "$ENV" $C Rscript 20b_contrast_screen.Rscript
cd $D/auk_ms/figure_code
apptainer exec --env "$ENV" $C Rscript fig3_genomewide.R
