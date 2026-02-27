#!/bin/bash
#SBATCH --job-name=fet_manhattan
#SBATCH --partition=128x24
#SBATCH --time=01:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=1
#SBATCH --output=fet.out
#SBATCH --error=fet.err

set -euo pipefail

# This job reads a PoPoolation2 FET file, computes cumulative genome positions using a .fai index
# (scaffolds ordered largest → smallest), applies BH-FDR, and writes a Manhattan-style plot to a PNG.

echo "[$(date)] Starting FET Manhattan plot job on $(hostname)"

source /hb/home/aschjeld/miniconda3/etc/profile.d/conda.sh
conda activate /hb/home/aschjeld/miniconda3/envs/R_env

Rscript -e "library(data.table); library(ggplot2); cat('R packages OK\n')"
Rscript FETman.r

echo "[$(date)] Finished. Output PNG should be in the path defined inside FETman.r"
