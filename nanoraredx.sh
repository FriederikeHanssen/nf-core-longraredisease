#!/bin/bash
#SBATCH --cpus-per-task=32
#SBATCH --time 4-0:0:0
#SBATCH --account beggsa-genomicsbirmingham
#SBATCH --qos bbshort
#SBATCH --time 10

module purge; module load bluebear
module load bear-apps/2022b
module load Nextflow/24.04.2

nextflow run main.nf -profile test,singularity