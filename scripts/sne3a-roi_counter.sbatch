#!/bin/bash
#SBATCH -p normal
#SBATCH -c 1
#SBATCH --cpus-per-task=4
#SBATCH --time 00:30:00
#SBATCH --output=/home/groups/rairan/NODDI/slurm_outputs//slurm-%j.out
#

rootdir=${1}
d=${2}
scripts_dir=${3}
miracl_container_path=${4}
lut_path=${5}
labels_path=${6}
output_dir=${7}
organism=${8}

cd $rootdir

mask_img=$(find "${rootdir}/${d}" -type f -name "template_mask.nii.gz")

singularity exec $miracl_container_path fslmaths $mask_img -bin "${rootdir}/${d}/brain_mask_orig.nii.gz"
singularity exec $miracl_container_path fslmaths "${rootdir}/${d}/brain_mask_orig.nii.gz" -ero "${rootdir}/${d}/brain_mask.nii.gz"
rm "${rootdir}/${d}/brain_mask_orig.nii.gz"

singularity exec $miracl_container_path python ${scripts_dir}/sne3b-roi_counter.py $d $lut_path $labels_path $output_dir $organism

