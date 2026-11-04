#!/bin/bash
#SBATCH --job-name=trp6_pm6_100fs
#SBATCH --partition=32core_partition
#SBATCH --nodelist=node3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=qmmm_pm6_100fs.slurm-%j.out

set -euo pipefail

module load amber/26.0

sander -O \
  -i qmmm_pm6_100fs.in \
  -o qmmm_pm6_100fs.out \
  -p ../step3_input.parm7 \
  -c eq_0.rst7 \
  -r qmmm_pm6_100fs.rst7 \
  -x qmmm_pm6_100fs.nc \
  -inf qmmm_pm6_100fs.mdinfo
