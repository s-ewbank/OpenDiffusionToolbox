#!/bin/bash
#SBATCH -p normal
#SBATCH -c 1
#SBATCH --cpus-per-task=4
#SBATCH --time 06:00:00
#SBATCH --output=/home/groups/rairan/NODDI/slurm_outputs//slurm-%j.out
#

scripts_dir=${1}
miracl_container_path=${2}
miracl_output_dir=${3}
mask_path=${4} 
baseline_tmap_dir=${5} 
avg_vol_dir=${6}
groups_file=${7}
organism=${8} 
opt1=${9}
opt2=${10}
opt3=${11}
volume_types=("${@:12}")

for vol in "${volume_types[@]}"
do
    mkdir -p ${baseline_tmap_dir}/tmp/
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Making t-maps for ${vol}"
    singularity exec $miracl_container_path python ${scripts_dir}/sne4c-tmaps.py ${miracl_output_dir} ${mask_path} ${baseline_tmap_dir} ${avg_vol_dir} ${groups_file} ${organism} ${opt1} ${opt2} ${vol} ${opt3}
    rm -r ${baseline_tmap_dir}/tmp/
done