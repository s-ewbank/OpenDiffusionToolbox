#!/bin/bash
#SBATCH -p normal
#SBATCH -c 1
#SBATCH --cpus-per-task=4
#SBATCH --time 4:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=snewbank@stanford.edu
#SBATCH --output=/home/groups/rairan/NODDI/slurm_outputs/slurm-%j.out
#

rootdir=${1}
d=${2}
scripts_dir=${3}
mdt_container_path=${4}
miracl_container_path=${5}

#######################################
# SET VARIABLES
#######################################
protocol="protocol.prtcl"

cd $rootdir
ml biology
ml fsl
ml system
ml clinfo


if [ ! -d "batch_output" ]; then mkdir batch_output; fi

#if [ ! -f "${rootdir}/filt_log.txt" ]; then touch "${rootdir}/filt_log.txt"; fi

subdirs=$(ls)
n_subdirs=`find . -mindepth 1 -maxdepth 1 -type d | wc -l`
    
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Variables set as:"
echo " -- Input data:"
echo "    > Root directory: $rootdir"
echo "    > Working directory: $d"
echo "    > Number subdirectories with presumed DWI volumes: $n_subdirs"
echo "    > List of subdirectories:"
for subdir in $subdirs; do echo "      + $subdir"; done
    
#################################################################################
# EDDY CORRECT, MASK GENERATION, THEN DTIFIT
#################################################################################
cd $d

singularity exec $miracl_container_path python ${scripts_dir}/sne1b-filt-mppca.py ${rootdir}/${d}


if [ ! -d "output" ]; then mkdir output; fi


for slice in slices/filt_*.nii.gz; do

    if [ ! -d "output/dti_output_dir" ]; then mkdir "output/dti_output_dir"; fi
    if [ ! -d "output/noddi_output_dir" ]; then mkdir "output/noddi_output_dir"; fi
    
    slice_i="${slice##*filt_}"
    slice_i="${slice_i%%.nii.gz*}"

    mask=slices/mask_${slice_i}.nii.gz
    dwi=slices/filt_${slice_i}.nii.gz
    bval=slices/bval_${slice_i}.bval
    bvec=slices/bvec_${slice_i}.bvec
    protocol=slices/protocol_${slice_i}.prtcl
    
    #######################################
    # DTIFIT
    #######################################
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Fitting DTI to data with dtifit and outputting FA, MD, etc. in directory '${rootdir}/output/dti_output_dir/' with dtifit..."
    echo ""
    dtifit --data=$dwi --mask=$mask --bvecs=$bvec --bvals=$bval --kurt --kurtdir --out=output/dti_output_dir/dti
    
    #######################################
    # RUN MDT
    #######################################
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Done creating mask '$mask' and fitting DTI, now creating protocol $protocol in subdirectory '$d'..."
    echo ""
    singularity exec $mdt_container_path env PYOPENCL_NO_CACHE=1 mdt-create-protocol $bvec $bval -o $protocol
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Done creating '$protocol', now fitting NODDI to DWI volume $dwi..."
    singularity exec $mdt_container_path env PYOPENCL_NO_CACHE=1 mdt-model-fit -o output/noddi_output_dir/ "NODDI" $dwi $protocol $mask
    echo ""
    
    if [ ! -d "output/final_output_dir" ]; then 
    mkdir "output/final_output_dir"
    cp output/dti_output_dir/dti_FA.nii.gz output/final_output_dir/FA.nii.gz
    cp output/dti_output_dir/dti_MD.nii.gz output/final_output_dir/MD.nii.gz
    cp output/dti_output_dir/dti_S0.nii.gz output/final_output_dir/S0.nii.gz
    cp output/noddi_output_dir/NODDI/ODI.nii.gz output/final_output_dir/ODI.nii.gz
    cp output/noddi_output_dir/NODDI/NDI.nii.gz output/final_output_dir/NDI.nii.gz
    rm -r output/noddi_output_dir
    rm -r output/dti_output_dir
    else
    fslmaths output/final_output_dir/FA.nii.gz -add output/dti_output_dir/dti_FA.nii.gz output/final_output_dir/FA.nii.gz
    fslmaths output/final_output_dir/MD.nii.gz -add output/dti_output_dir/dti_MD.nii.gz output/final_output_dir/MD.nii.gz
    fslmaths output/final_output_dir/S0.nii.gz -add output/dti_output_dir/dti_S0.nii.gz output/final_output_dir/S0.nii.gz
    fslmaths output/final_output_dir/ODI.nii.gz -add output/noddi_output_dir/NODDI/ODI.nii.gz output/final_output_dir/ODI.nii.gz
    fslmaths output/final_output_dir/NDI.nii.gz -add output/noddi_output_dir/NODDI/NDI.nii.gz output/final_output_dir/NDI.nii.gz
    rm -r output/noddi_output_dir
    rm -r output/dti_output_dir
    fi
    
done

mv slices/filt_log.csv output/final_output_dir/filt_log.csv
mv output/final_output_dir ${d}_output
mv ${d}_output ../batch_output/${d}_output
rm -r slices






