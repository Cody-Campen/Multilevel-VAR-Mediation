#!/bin/bash
#SBATCH --qos=normal
#SBATCH --partition=basic
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8gb
#SBATCH --time=2-00:00:00
#SBATCH --job-name=informative=FALSE_diffuseness=1
#SBATCH --error=%u_%x.out
#SBATCH --output=%u_%x.out
SYSTEM=$SLURM_SUBMIT_DIR/System
cd /storage/work/cac7655/Multilevel-VAR-Mediation
module load r/4.5.0
module load lapack/3.11.0-gcc
export PKG_CONFIG_PATH=/storage/work/cac7655/JAGS-4.3.2/etc
export LD_LIBRARY_PATH=/storage/work/cac7655/JAGS-4.3.2/lib:$LD_LIBRARY_PATH
export R_LIBS=/storage/work/cac7655/RLibs/4.5.0
R <  runsim_seed=27_informative=FALSE_diffuseness=1.R  --no-save
