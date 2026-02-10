#!/bin/bash
#SBATCH -p normal
#SBATCH -c 1
#SBATCH --cpus-per-task=4
#SBATCH --time 3:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=snewbank@stanford.edu
#SBATCH --output=/scratch/groups/rairan/NODDI/test/slurm-%j.out

# Currently diffusiontoolkit-DEV.sif is based on Singularity12.def
# And diffusiontoolkit.sif is based on Singularity11.def
# Singularity_FINAL-v2 edit is to accept conda terms of service (otherwise get nonresponse error)

singularity build /scratch/groups/rairan/NODDI/test/diffusiontoolkit-FINAL.sif /home/groups/rairan/OpenDiffusionToolbox/singularity_recipe/Singularity_FINAL-v2.def
