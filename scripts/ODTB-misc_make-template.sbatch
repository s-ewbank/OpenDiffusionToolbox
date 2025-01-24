#!/bin/bash
#SBATCH -p normal
#SBATCH -c 1
#SBATCH --cpus-per-task=4
#SBATCH --time 48:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=snewbank@stanford.edu
#SBATCH --output=/scratch/users/snewbank/slurm_output//slurm-%j.out
#

script_dir_name=${1}
rootdir=${2}
vol_dir=${3}
miracl_container_path=${4}
step=${5}
volumes=${6}

echo "[$(date '+%Y-%m-%d %H:%M:%S')] RUNNING REGISTRATION SCRIPT"
echo $script_dir_name
echo " "

ml biology
ml fsl
ml ants

export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$SLURM_CPUS_PER_TASK

cd $rootdir
cd $vol_dir
mkdir ${vol_dir}_reg_output
 
#######################################
# MIRACL - SORTING OUT VOLUMES
#######################################
mkdir mr_allen_reg

cp FA.nii.gz ${vol_dir}_reg_output/${vol_dir}_FA.nii.gz
cp MD.nii.gz ${vol_dir}_reg_output/${vol_dir}_MD.nii.gz
cp NDI.nii.gz ${vol_dir}_reg_output/${vol_dir}_NDI.nii.gz
cp ODI.nii.gz ${vol_dir}_reg_output/${vol_dir}_ODI.nii.gz
t2_dat=S0.nii.gz
cp $t2_dat ${vol_dir}_reg_output/${vol_dir}_S0.nii.gz

FA_vol=`find ${vol_dir}_reg_output/*_FA.nii.gz`
MD_vol=`find ${vol_dir}_reg_output/*_MD.nii.gz`
ODI_vol=`find ${vol_dir}_reg_output/*_ODI.nii.gz`
NDI_vol=`find ${vol_dir}_reg_output/*_NDI.nii.gz`
S0_vol=`find ${vol_dir}_reg_output/*_S0.nii.gz`
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting MIRACL registration for ${vol_dir}."
echo "Volumes to be registered are:"
echo "   - FA volume: $FA_vol"
echo "   - MD volume: $MD_vol"
echo "   - NDI volume: $NDI_vol"
echo "   - ODI volume: $ODI_vol"
echo " "

## IF on step 1

if [[ "$step" == 1 ]]; then
    #######################################
    # MIRACL - DOING ANTS ON T2
    #######################################
        
    mkdir reg_final
    
    pre_reg_mask=pre_reg_mask.nii.gz
    
    mkdir ${rootdir}/${vol_dir}/temp
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Bias correcting"
    # FAST for bias correcting 
    fast -l 1 -I 20 -B $t2_dat
    rm ${t2_dat//.nii.gz}_seg.nii.gz ${t2_dat//.nii.gz}_pve*.nii.gz ${t2_dat//.nii.gz}_mixeltype.nii.gz
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Extracting brain"
    # Take square root and do initial thresholding at 40th percentile
    fslmaths ${t2_dat//.nii.gz}_restore.nii.gz -sqrt ${rootdir}/${vol_dir}/temp/sqrt.nii.gz
    fslmaths ${rootdir}/${vol_dir}/temp/sqrt.nii.gz -thrP 40 ${rootdir}/${vol_dir}/temp/thresh.nii.gz
    
    # Binarize, smooth, threshold and binarize to get a mask excluding peripheral voxels which have survived thresholding, then smooth to create a cloud and multiply by original image
    fslmaths ${rootdir}/${vol_dir}/temp/thresh.nii.gz -bin -s 1.5 -thrP 80 -bin -s 1.5 ${rootdir}/${vol_dir}/temp/cloud.nii.gz
    fslmaths ${rootdir}/${vol_dir}/temp/cloud.nii.gz -mul ${t2_dat//.nii.gz}_restore.nii.gz -thrP 20 ${rootdir}/${vol_dir}/temp/cloud_intermed.nii.gz
    fslmaths ${rootdir}/${vol_dir}/temp/cloud_intermed.nii.gz -bin -ero -ero -ero -s 0.5 -thr 0.1 -bin -dilD -dilD ${rootdir}/${vol_dir}/temp/cloud_mask.nii.gz
    fslmaths ${t2_dat//.nii.gz}_restore.nii.gz -mas ${rootdir}/${vol_dir}/temp/cloud_mask.nii.gz -thrP 30 -bin -ero -ero -dilD -dilD $pre_reg_mask
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Masking non-brain tissue and ventricles"
    fslmaths $MD_vol -mas $pre_reg_mask -thrP 75 -bin -mul -1 -add 1 ${rootdir}/${vol_dir}/ventricles.nii.gz
    fslmaths ${t2_dat//.nii.gz}_restore.nii.gz -mas $pre_reg_mask -mas ${rootdir}/${vol_dir}/ventricles.nii.gz ${t2_dat//.nii.gz}_masked.nii.gz
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Clip at 90th percentile and apply POW!!"
    p90=$(fslstats ${t2_dat//.nii.gz}_masked.nii.gz -P 90 | bc)  
    fslmaths ${t2_dat//.nii.gz}_masked.nii.gz -min $p90 -pow 2 ${t2_dat//.nii.gz}_masked_min_pow.nii.gz
    
    rm -r ${rootdir}/${vol_dir}/temp
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Done preprocessing S0 image - starting FLIRT reg step"
        
    singularity exec $miracl_container_path cp /code/atlases/ara/template/average_template_50um.nii.gz ${rootdir}/${vol_dir}/template.nii.gz
    fslmaths ${rootdir}/${vol_dir}/template.nii.gz -subsamp2 ${rootdir}/${vol_dir}/template.nii.gz
    fslmaths ${rootdir}/${vol_dir}/template.nii.gz -thrP 5 -bin -dilD -dilD -dilD -dilD ${rootdir}/${vol_dir}/template_mask.nii.gz
    
    flirt -in ${t2_dat//.nii.gz}_masked_min_pow.nii.gz \
        -ref ${rootdir}/${vol_dir}/template.nii.gz \
        -init ${script_dir_name}/mouse_S0_to_ref.mat \
        -searchrx -15 15 -searchry -15 15 -searchrz -15 15 \
        -omat ${rootdir}/${vol_dir}/subj_S0_to_ref.mat \
        -noresampblur \
        -out ${t2_dat//.nii.gz}_flirt.nii.gz
    
    fslmaths ${t2_dat//.nii.gz}_flirt.nii.gz -mas ${rootdir}/${vol_dir}/template_mask.nii.gz ${t2_dat//.nii.gz}_flirt.nii.gz
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Done with flirt reg step - starting ANTs reg step"
    
    antsRegistration --dimensionality 3 --float 0 \
        --output [${rootdir}/${vol_dir}/mr_allen_reg/allen_mr_ants,${rootdir}/${vol_dir}/mr_allen_reg/allen_mr_antsWarped.nii.gz,${rootdir}/${vol_dir}/mr_allen_reg/allen_mr_antsInverseWarped.nii.gz] \
        --interpolation Linear --use-histogram-matching 0 --winsorize-image-intensities [0.005,0.995] \
        --initial-moving-transform [${t2_dat//.nii.gz}_flirt.nii.gz,${rootdir}/${vol_dir}/template.nii.gz,1] \
        --transform Affine[0.1] \
        --metric CC[${t2_dat//.nii.gz}_flirt.nii.gz,${rootdir}/${vol_dir}/template.nii.gz,1,8] \
        --convergence [2000x700x550x200,1e-10,15] --shrink-factors 4x2x1x1 --smoothing-sigmas 3x2x1x0vox \
        --transform BSplineSyN[0.1,26,0,3] \
        --metric CC[${t2_dat//.nii.gz}_flirt.nii.gz,${rootdir}/${vol_dir}/template.nii.gz,1,8] \
        --convergence [100x70x50x20,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox
        
    #antsRegistration --dimensionality 3 --float 0 \
    #    --output [${rootdir}/${vol_dir}/mr_allen_reg/allen_mr_ants,${rootdir}/${vol_dir}/mr_allen_reg/allen_mr_antsWarped.nii.gz,${rootdir}/${vol_dir}/mr_allen_reg/allen_mr_antsInverseWarped.nii.gz] \
    #    --interpolation Linear --use-histogram-matching 0 --winsorize-image-intensities [0.005,0.995] \
    #    --initial-moving-transform [${t2_dat//.nii.gz}_denoised_scaled_flirt_noventricles.nii.gz,${rootdir}/${vol_dir}/template.nii.gz,1] \
    #    --transform Rigid[0.1] \
    #    --metric MI[${t2_dat//.nii.gz}_masked_min_pow.nii.gz,${rootdir}/${vol_dir}/template.nii.gz,1,2] \
    #    --convergence [2000x700x550x200,1e-10,15] --shrink-factors 4x2x1x1 --smoothing-sigmas 3x2x1x0vox --transform Affine[0.1] \
    #    --metric CC[${t2_dat//.nii.gz}_masked_min_pow.nii.gz,${rootdir}/${vol_dir}/template.nii.gz,1,8] \
    #    --convergence [2000x700x550x200,1e-10,15] --shrink-factors 4x2x1x1 --smoothing-sigmas 3x2x1x0vox --transform BSplineSyN[0.1,26,0,3] \
    #    --metric CC[${t2_dat//.nii.gz}_masked_min_pow.nii.gz,${rootdir}/${vol_dir}/template.nii.gz,1,8] \
    #    --convergence [100x70x50x20,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox
        
        
    
    #######################################
    # MIRACL - APPLYING ANTS TO OTHER VOLS
    #######################################
    
    volumes=("$FA_vol" "$MD_vol" "$ODI_vol" "$NDI_vol" "$S0_vol")
    
    for vol in "${volumes[@]}"
    do
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Working on volume ${vol}"
        
        flirt -in $vol \
            -ref ${rootdir}/${vol_dir}/template.nii.gz \
            -init ${rootdir}/${vol_dir}/subj_S0_to_ref.mat \
            -applyxfm -noresampblur \
            -out ${vol//.nii.gz}_flirt.nii.gz
        
        flirt -in ${rootdir}/${vol_dir}/mr_allen_reg/allen_mr_ants1InverseWarp.nii.gz -ref ${rootdir}/${vol_dir}/template.nii.gz -applyxfm -out ${rootdir}/${vol_dir}/mr_allen_reg/allen_mr_ants1InverseWarp_templateres.nii.gz
        
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Registering"
        singularity exec $miracl_container_path antsApplyTransforms -d 3 \
            -r ${rootdir}/${vol_dir}/template.nii.gz \
            -i ${vol//.nii.gz}_flirt.nii.gz -n Linear \
            -t [${rootdir}/${vol_dir}/mr_allen_reg/allen_mr_ants0GenericAffine.mat,1] [${rootdir}/${vol_dir}/mr_allen_reg/allen_mr_ants1InverseWarp.nii.gz,0] \
            -o ${vol%.nii.gz}_reg.nii.gz --float -v
        
        #fslmaths ${vol%.nii.gz}_reg.nii.gz -subsamp2 -subsamp2 ${vol%.nii.gz}_reg_ss2.nii.gz
        
    done
    
    mv ${vol_dir}_reg_output ../batch_reg_output
    cd ..

fi






