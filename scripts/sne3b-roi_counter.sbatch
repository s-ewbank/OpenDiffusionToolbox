#!/bin/bash
#SBATCH -p normal
#SBATCH -c 1
#SBATCH --cpus-per-task=4
#SBATCH --time 00:30:00
#SBATCH --output=/home/groups/rairan/NODDI/slurm_outputs/slurm-%j.out
#

rootdir=${1}
d=${2}
lut_path=${3}
labels_path=${4}
output_dir=${5}

cd $rootdir
sing_path=/scratch/groups/rairan/NODDI/MIRACL/miracl_latest.sif

ndi_img=$(find "${rootdir}/${d}" -type f -name "*_NDI_reg.nii.gz")

singularity exec $sing_path fslmaths $ndi_img -bin "${rootdir}/${d}/miracl_brain_mask_orig.nii.gz"
singularity exec $sing_path fslmaths "${rootdir}/${d}/miracl_brain_mask_orig.nii.gz" -ero "${rootdir}/${d}/miracl_brain_mask.nii.gz"
rm "${rootdir}/${d}/miracl_brain_mask_orig.nii.gz"

singularity exec $sing_path python /home/groups/rairan/NODDI/sne3d-roi_counter.py $d $lut_path $labels_path $output_dir

