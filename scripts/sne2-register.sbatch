#!/bin/bash
#SBATCH -p normal
#SBATCH -c 1
#SBATCH --cpus-per-task=4
#SBATCH --time 4:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=snewbank@stanford.edu
#SBATCH --output=/home/groups/rairan/NODDI/slurm_outputs//slurm-%j.out
#

script_dir_name=${1}
rootdir=${2}
vol_dir=${3}
do_ants=${4}
apply_ants_to_all_vols=${5}
miracl_container_path=${6}
organism=${7}
volumes=${8}

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

if [[ "$organism" == "mouse" ]]; then
    #######################################
    # SET VARIABLES AND LOAD MODULES
    #######################################
    labels=/code/atlases/ara/annotation/annotation_hemi_combined_50um.nii.gz
     
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
    
    #######################################
    # MIRACL - DOING ANTS ON T2
    #######################################
    
    if [[ $do_ants -ne 0 ]]
    then
        
        mkdir reg_final
    
        mkdir tmp
        num_slices=$(fslhd $t2_dat | awk '/^dim2/{print $2}')
        half_slices=$((num_slices / 2))
        
        fslroi $t2_dat tmp/bottom_half.nii.gz 0 -1 0 $half_slices 0 -1
        fslmaths tmp/bottom_half.nii.gz -sub 100000000 -bin tmp/bottom_half.nii.gz
        
        fslroi $t2_dat tmp/top_half.nii.gz 0 -1 $half_slices -1 0 -1
        fslmaths tmp/top_half.nii.gz -add 1000 -bin tmp/top_half.nii.gz
        fslmerge -y tmp/combined_image.nii.gz tmp/bottom_half.nii.gz tmp/top_half.nii.gz
        fslmaths tmp/combined_image.nii.gz -s 3 -mul 5 -mul $t2_dat -add $t2_dat ${t2_dat//.nii.gz}_grad.nii.gz
        
        rm -r tmp
        
        bash ${script_dir_name}/sne2b-rbet.sh $rootdir $vol_dir ${t2_dat//.nii.gz}_grad.nii.gz ${rootdir}/${vol_dir}/pre_reg_mask.nii.gz
        fslmaths $t2_dat -mas ${rootdir}/${vol_dir}/pre_reg_mask.nii.gz -s 0.1 ${t2_dat//.nii.gz}_masked.nii.gz
        
        #get range of original image
        min=$(fslstats ${t2_dat//.nii.gz}_grad_masked.nii.gz -P 0)
        max=$(fslstats ${t2_dat//.nii.gz}_grad_masked.nii.gz -P 100)
        sc=$(echo "$max - $min" | bc)
        
        #Extract high frequency noise-type voxels and replace with smoothed voxels
        fslmaths ${t2_dat//.nii.gz}_grad_masked.nii.gz -s 0.2 ${t2_dat//.nii.gz}_smoothed.nii.gz
        fslmaths ${t2_dat//.nii.gz}_smoothed.nii.gz -sub ${t2_dat//.nii.gz}_grad_masked.nii.gz ${t2_dat//.nii.gz}_diff.nii.gz
        lo=$(fslstats ${t2_dat//.nii.gz}_diff.nii.gz -P 10)
        hi=$(fslstats ${t2_dat//.nii.gz}_diff.nii.gz -P 90)
        fslmaths ${t2_dat//.nii.gz}_diff.nii.gz -thr $lo -uthr $hi -bin ${t2_dat//.nii.gz}_diff_mask.nii.gz
        fslmaths ${t2_dat//.nii.gz}_diff_mask.nii.gz -s 0.2 -uthr 0.6 -bin -mas ${t2_dat//.nii.gz}_diff_mask.nii.gz ${t2_dat//.nii.gz}_diff_sub_mask.nii.gz
        
        #also do nonlinear scaling
        fslmaths ${t2_dat//.nii.gz}_diff.nii.gz -mas ${t2_dat//.nii.gz}_diff_sub_mask.nii.gz -mul 0.7 -add ${t2_dat//.nii.gz}_grad_masked.nii.gz ${t2_dat//.nii.gz}_denoised.nii.gz
        min2=$(fslstats ${t2_dat//.nii.gz}_denoised.nii.gz -P 0)
        max2=$(fslstats ${t2_dat//.nii.gz}_denoised.nii.gz -P 100)
        sc2=$(echo "$max2 - $min2" | bc)
        
        #rescale to original range
        fslmaths ${t2_dat//.nii.gz}_denoised.nii.gz -sub $min2 -div $sc2 -mul $sc -add $min ${t2_dat//.nii.gz}_denoised_scaled.nii.gz
        
        rm ${t2_dat//.nii.gz}_denoised.nii.gz ${t2_dat//.nii.gz}_diff_sub_mask.nii.gz ${t2_dat//.nii.gz}_diff_mask.nii.gz ${t2_dat//.nii.gz}_diff.nii.gz ${t2_dat//.nii.gz}_smoothed.nii.gz
        
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running mri_allen_ants"
        #singularity exec $miracl_container_path miracl reg mri_allen_ants -i ${t2_dat//.nii.gz}_denoised_scaled.nii.gz \
        #    -o RSP -m combined -v 50 -b 1 -s 0
            
        singularity exec $miracl_container_path cp /code/atlases/ara/template/average_template_50um_OBmasked.nii.gz ${rootdir}/${vol_dir}/template.nii.gz
        fslmaths ${rootdir}/${vol_dir}/template.nii.gz -subsamp2 -subsamp2 ${rootdir}/${vol_dir}/template_ss4.nii.gz
        
        flirt -in ${t2_dat//.nii.gz}_denoised_scaled.nii.gz \
            -ref ${rootdir}/${vol_dir}/template_ss4.nii.gz \
            -init ${script_dir_name}/mouse_S0_to_ref.mat \
            -searchrx -15 15 -searchry -15 15 -searchrz -15 15 \
            -omat ${rootdir}/${vol_dir}/subj_S0_to_ref.mat \
            -noresampblur \
            -out ${t2_dat//.nii.gz}_denoised_scaled_flirt.nii.gz
        
        flirt -in $MD_vol \
            -ref ${rootdir}/${vol_dir}/template_ss4.nii.gz \
            -init ${rootdir}/${vol_dir}/subj_S0_to_ref.mat \
            -applyxfm -noresampblur \
            -out ${MD_vol//.nii.gz}_flirt.nii.gz
            
        fslmaths ${MD_vol//.nii.gz}_flirt.nii.gz -thrP 80 -bin -mul -1 -add 1 ${rootdir}/${vol_dir}/ventricles.nii.gz
        min=$(fslstats ${t2_dat//.nii.gz}_denoised_scaled_flirt.nii.gz -P 1)
        min_float=$(echo $min | bc)
        fslmaths ${t2_dat//.nii.gz}_denoised_scaled_flirt.nii.gz -sub $min_float -mas ${rootdir}/${vol_dir}/ventricles.nii.gz -pow 1.5 ${t2_dat//.nii.gz}_denoised_scaled_flirt_noventricles.nii.gz
                
        #fslmaths ${rootdir}/${vol_dir}/mr_allen_reg/allen_mr_antsInverseWarped.nii.gz -thr 0 -pow 3 -subsamp2 -subsamp2 ${rootdir}/${vol_dir}/mr_allen_reg/round2_input.nii.gz
        #
        antsRegistration --dimensionality 3 --float 0 \
            --output [${rootdir}/${vol_dir}/mr_allen_reg/allen_mr_ants,${rootdir}/${vol_dir}/mr_allen_reg/allen_mr_antsWarped.nii.gz,${rootdir}/${vol_dir}/mr_allen_reg/allen_mr_antsInverseWarped.nii.gz] \
            --interpolation Linear --use-histogram-matching 0 --winsorize-image-intensities [0.005,0.995] \
            --initial-moving-transform [${t2_dat//.nii.gz}_denoised_scaled_flirt_noventricles.nii.gz,${rootdir}/${vol_dir}/template.nii.gz,1] \
            --transform Rigid[0.1] \
            --metric MI[${t2_dat//.nii.gz}_denoised_scaled_flirt_noventricles.nii.gz,${rootdir}/${vol_dir}/template.nii.gz,1,2] \
            --convergence [2000x700x550x200,1e-10,15] --shrink-factors 4x2x1x1 --smoothing-sigmas 3x2x1x0vox --transform Affine[0.1] \
            --metric CC[${t2_dat//.nii.gz}_denoised_scaled_flirt_noventricles.nii.gz,${rootdir}/${vol_dir}/template.nii.gz,1,8] \
            --convergence [2000x700x550x200,1e-10,15] --shrink-factors 4x2x1x1 --smoothing-sigmas 3x2x1x0vox --transform BSplineSyN[0.1,26,0,3] \
            --metric CC[${t2_dat//.nii.gz}_denoised_scaled_flirt_noventricles.nii.gz,${rootdir}/${vol_dir}/template.nii.gz,1,8] \
            --convergence [100x70x50x20,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox
            
    fi
    
    #######################################
    # MIRACL - APPLYING ANTS TO OTHER VOLS
    #######################################
    if [[ $apply_ants_to_all_vols -ne 0 ]]
    then
        #mask=${rootdir}/${vol_dir}/${vol_dir}_reg_output/miracl_mask.nii.gz
        preproc_1=${rootdir}/${vol_dir}/${vol_dir}_reg_output/preproc1.nii.gz
        preproc_2=${rootdir}/${vol_dir}/${vol_dir}_reg_output/preproc2.nii.gz
        preproc_3=${rootdir}/${vol_dir}/${vol_dir}_reg_output/preproc3.nii.gz
        preproc_4=${rootdir}/${vol_dir}/${vol_dir}_reg_output/preproc4.nii.gz
        preproc_4_masked=${rootdir}/${vol_dir}/${vol_dir}_reg_output/preproc4_masked.nii.gz
    
        #fslmaths ${rootdir}/${vol_dir}/mr_allen_reg/mr.nii.gz -bin $mask
        #fslmaths ${t2_dat//.nii.gz}_grad_masked.nii.gz -bin $mask
        
        volumes=("$FA_vol" "$MD_vol" "$ODI_vol" "$NDI_vol" "$S0_vol")
        
        for vol in "${volumes[@]}"
        do
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Working on volume ${vol}"
            
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Doing reorient steps"
            singularity exec $miracl_container_path c3d $vol -orient RSP -o $preproc_1
            singularity exec $miracl_container_path c3d $preproc_1 -spacing 1.500000x1.495330x4.000000mm -o $preproc_2
            singularity exec $miracl_container_path c3d $preproc_2 -spacing .150000x.149533x.400000mm -o $preproc_3
            singularity exec $miracl_container_path flirt -in $preproc_3 -ref ${rootdir}/${vol_dir}/mr_allen_reg/mr.nii.gz -out $preproc_4 -applyxfm -usesqform 
            
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Masking"
            #fslmaths $preproc_4 -mas $mask $preproc_4_masked
            
            flirt -in $vol \
                -ref ${rootdir}/${vol_dir}/template.nii.gz \
                -init ${rootdir}/${vol_dir}/subj_S0_to_ref.mat \
                -applyxfm -noresampblur \
                -out ${vol//.nii.gz}_flirt.nii.gz
            
            flirt -in ${rootdir}/${vol_dir}/mr_allen_reg/allen_mr_ants1InverseWarp.nii.gz -ref ${rootdir}/${vol_dir}/template.nii.gz -applyxfm -out ${rootdir}/${vol_dir}/mr_allen_reg/allen_mr_ants1InverseWarp_templateres.nii.gz
            
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Registering"
            #change i to $preproc_4 if needed
            singularity exec $miracl_container_path antsApplyTransforms -d 3 \
                -r ${rootdir}/${vol_dir}/template.nii.gz \
                -i ${vol//.nii.gz}_flirt.nii.gz -n Linear \
                -t [${rootdir}/${vol_dir}/mr_allen_reg/allen_mr_ants0GenericAffine.mat,1] [${rootdir}/${vol_dir}/mr_allen_reg/allen_mr_ants1InverseWarp_templateres.nii.gz,0] \
                -o ${vol%.nii.gz}_reg.nii.gz --float -v
            
            #fslmaths ${vol%.nii.gz}_reg.nii.gz -subsamp2 -subsamp2 ${vol%.nii.gz}_reg_ss2.nii.gz
            #
            #singularity exec $miracl_container_path antsApplyTransforms -d 3 \
            #    -r ${rootdir}/${vol_dir}/template.nii.gz \
            #    -i ${vol%.nii.gz}_reg_ss2.nii.gz -n Linear \
            #    -t [${rootdir}/${vol_dir}/mr_allen_reg/round2_allen_mr_ants0GenericAffine.mat,1] [${rootdir}/${vol_dir}/mr_allen_reg/round2_allen_mr_ants1InverseWarp.nii.gz,0] \
            #    -o ${vol%.nii.gz}_reg2.nii.gz --float -v
                
            #flirt -in ${rootdir}/${vol_dir}/mr_allen_reg/round2_allen_mr_ants1InverseWarp.nii.gz -ref ${rootdir}/${vol_dir}/template.nii.gz \
            #    -applyxfm -out ${rootdir}/${vol_dir}/mr_allen_reg/round2_allen_mr_ants1InverseWarp_us4.nii.gz
            #
            #singularity exec $miracl_container_path antsApplyTransforms -d 3 \
            #    -r ${rootdir}/${vol_dir}/template.nii.gz \
            #    -i ${vol%.nii.gz}_reg.nii.gz -n Linear \
            #    -t [${rootdir}/${vol_dir}/mr_allen_reg/round2_allen_mr_ants0GenericAffine.mat,1] [${rootdir}/${vol_dir}/mr_allen_reg/round2_allen_mr_ants1InverseWarp_us4.nii.gz,0] \
            #    -o ${vol%.nii.gz}_reg3.nii.gz --float -v
                
            rm $preproc_1
            rm $preproc_2
            rm $preproc_3
            rm $preproc_4
            #mv $preproc_4_masked ${vol%.nii.gz}_ss.nii.gz
            
        done
    
    fi
    
fi


if [[ "$organism" == "human" ]]; then
    
    #######################################
    # SET VARIABLES AND LOAD MODULES
    #######################################
    
    ref=$FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz
    cort_labels=$FSLDIR/data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr25-2mm.nii.gz
    subcort_labels=$FSLDIR/data/atlases/HarvardOxford/HarvardOxford-sub-maxprob-thr25-2mm.nii.gz
     
    #######################################
    # REG - SORTING OUT VOLUMES
    #######################################
    mkdir ${vol_dir}_reg_output
    
    cp FA.nii.gz ${vol_dir}_reg_output/${vol_dir}_FA.nii.gz
    cp MD.nii.gz ${vol_dir}_reg_output/${vol_dir}_MD.nii.gz
    cp NDI.nii.gz ${vol_dir}_reg_output/${vol_dir}_NDI.nii.gz
    cp ODI.nii.gz ${vol_dir}_reg_output/${vol_dir}_ODI.nii.gz
    cp S0.nii.gz ${vol_dir}_reg_output/${vol_dir}_S0.nii.gz
    
    FA_vol=`find ${vol_dir}_reg_output/*_FA.nii.gz`
    MD_vol=`find ${vol_dir}_reg_output/*_MD.nii.gz`
    ODI_vol=`find ${vol_dir}_reg_output/*_ODI.nii.gz`
    NDI_vol=`find ${vol_dir}_reg_output/*_NDI.nii.gz`
    S0_vol=`find ${vol_dir}_reg_output/*_S0.nii.gz`
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting registration for ${vol_dir}."
    echo "Volumes to be registered are:"
    echo "   - FA volume: $FA_vol"
    echo "   - MD volume: $MD_vol"
    echo "   - NDI volume: $NDI_vol"
    echo "   - ODI volume: $ODI_vol"
    echo " "
    
    #######################################
    # REG - DOING ANTS ON T2
    #######################################
    
    if [[ $do_ants -ne 0 ]]
    then
        
        fslmaths $S0_vol -thrP 10 ${S0_vol//.nii.gz}_thr.nii.gz
        bet ${S0_vol//.nii.gz}_thr.nii.gz ${S0_vol//.nii.gz}_bet.nii.gz -m
    
        mv ${S0_vol//.nii.gz}_bet_mask.nii.gz ${vol_dir}_reg_output/mask.nii.gz
        
    #    antsRegistration --dimensionality 3 \
    #        --output [${S0_vol//.nii.gz}_bet_ants,${S0_vol//.nii.gz}_bet_antsWarped.nii.gz,${S0_vol//.nii.gz}_bet_antsInverseWarped.nii.gz] \
    #        --interpolation Linear \
    #        --winsorize-image-intensities [0.005,0.995] \
    #        --use-histogram-matching 1 \
    #        --transform Affine[0.1] \
    #        --metric CC[$ref,${S0_vol//.nii.gz}_bet.nii.gz,1,4] \
    #        --convergence [1000x350x250x100,1e-10,15] --shrink-factors 4x2x1x1 --smoothing-sigmas 3x2x1x0vox \
    #        --transform SyN[0.5,10,0] \
    #        --metric CC[$ref,${S0_vol//.nii.gz}_bet.nii.gz,1,4] \
    #        --convergence [200x150x100x50,1e-6,25] \
    #        --shrink-factors 2x1x1x1 \
    #        --smoothing-sigmas 2x1x0.5x0.5vox
    
        flirt -in ${S0_vol//.nii.gz}_bet.nii.gz -ref $ref -omat ${vol_dir}_reg_output/s0_to_ref.mat -out ${S0_vol//.nii.gz}_bet_reg.nii.gz -dof 12
            
    fi
    
    #######################################
    # REG - APPLYING ANTS TO OTHER VOLS
    #######################################
    if [[ $apply_ants_to_all_vols -ne 0 ]]
    then
        
        fslmaths $ref -bin ${vol_dir}_reg_output/template_mask.nii.gz
    #    aff=${S0_vol//.nii.gz}_bet_ants0GenericAffine.mat
    #    warp=${S0_vol//.nii.gz}_bet_ants1Warp.nii.gz
        
        volumes=("$FA_vol" "$MD_vol" "$ODI_vol" "$NDI_vol" "$S0_vol")
        
        for vol in "${volumes[@]}"
        do
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Working on volume ${vol}"
        
            flirt -in $vol -ref $ref -init ${vol_dir}_reg_output/s0_to_ref.mat -applyxfm -out ${vol%.nii.gz}_reg.nii.gz
            
    #        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Applying registration"
    #        echo "  Using transforms:"
    #        echo "  - ${aff}"
    #        echo "  - ${warp}"
            
    #        antsApplyTransforms -d 3 \
    #            -i $vol \
    #            -r $ref \
    #            -t $aff \
    #            -t $warp \
    #            -o ${vol//.nii.gz}_reg.nii.gz
                
            fslmaths ${vol//.nii.gz}_reg.nii.gz -mas ${vol_dir}_reg_output/template_mask.nii.gz ${vol//.nii.gz}_reg_masked.nii.gz
            
        done
    
    fi
    
fi

if [[ "$organism" == "rat" ]]; then
    rootdir=${1}
    vol_dir=${2}
    do_ants=${3}
    apply_ants_to_all_vols=${4}
    make_split_imgs=${5}
    
    #######################################
    # SET VARIABLES AND LOAD MODULES
    #######################################
    T2_template=/scratch/groups/rairan/PH_StrokeSC_MRI/WHS_SD_rat_atlas_v4_pack/WHS_SD_rat_T2star_v1.01_clipped_ss4_masked.nii.gz
    init=/home/groups/rairan/PH_StrokeSC_MRI/NEW_AND_IMPROVED/iso_to_ref_CORRECT.mat
    
    #######################################
    # SORTING OUT VOLUMES
    #######################################
    
    cp *_ISO.nii.gz ${vol_dir}_reg_output/${vol_dir}_ISO.nii.gz
    cp *_FA.nii.gz ${vol_dir}_reg_output/${vol_dir}_FA.nii.gz
    cp *_MD.nii.gz ${vol_dir}_reg_output/${vol_dir}_MD.nii.gz
    cp *_AD.nii.gz ${vol_dir}_reg_output/${vol_dir}_AD.nii.gz
    cp *_RD.nii.gz ${vol_dir}_reg_output/${vol_dir}_RD.nii.gz
    cp *_QA.nii.gz ${vol_dir}_reg_output/${vol_dir}_QA.nii.gz
    
    T2_struct_vol=`find ${vol_dir}_reg_output/${vol_dir}_ISO.nii.gz`
    FA_vol=`find ${vol_dir}_reg_output/*_FA.nii.gz`
    MD_vol=`find ${vol_dir}_reg_output/*_MD.nii.gz`
    AD_vol=`find ${vol_dir}_reg_output/*_AD.nii.gz`
    RD_vol=`find ${vol_dir}_reg_output/*_RD.nii.gz`
    QA_vol=`find ${vol_dir}_reg_output/*_QA.nii.gz`
    
    pre_reg_mask=${rootdir}/${vol_dir}/${vol_dir}_reg_output/pre_reg_mask.nii.gz
    post_reg_mask=${rootdir}/${vol_dir}/${vol_dir}_reg_output/T2_template_mask.nii.gz
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting registration for ${vol_dir}."
    echo "Volumes to be registered are:"
    echo "   - FA volume: $FA_vol"
    echo "   - QA volume: $QA_vol"
    echo "   - MD volume: $MD_vol"
    echo "   - RD volume: $RD_vol"
    echo "   - AD volume: $AD_vol"
    echo "   - ISO volume (struct): $T2_struct_vol"
    echo " "
    
    #######################################
    # EXTRACT BRAIN
    #######################################
    
    bash /home/groups/rairan/PH_StrokeSC_MRI/NEW_AND_IMPROVED/sne2c-bet.sh $rootdir $vol_dir $T2_struct_vol $pre_reg_mask
    
    #######################################
    # MIRACL - APPLYING ANTS TO OTHER VOLS
    #######################################
    
    volumes=("$T2_struct_vol" "$FA_vol" "$AD_vol" "$RD_vol" "$MD_vol" "$QA_vol")
    fslmaths $T2_template -bin $post_reg_mask
    
    for vol in "${volumes[@]}"
    do
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Working on volume ${vol}"
        
        if [ "$vol" == "$T2_struct_vol" ]; then
            
            fslmaths $vol -mas $pre_reg_mask ${vol%.nii.gz}_ss.nii.gz
            
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Registering to template with FLIRT [for T2w only]"
            flirt -in ${vol%.nii.gz}_ss.nii.gz -ref $T2_template -omat ${vol_dir}_reg_output/iso_to_template.mat -init $init -out ${vol%.nii.gz}_flirt.nii.gz -dof 12 \
                -searchrx -5 5 -searchry -5 5 -searchrz -5 5
                
                
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Using ANTS for nonlinear reg - hopefully this does not take 5-ever [for T2w only]"
            antsRegistration --dimensionality 3 \
                --output [${rootdir}/${vol_dir}/${vol//.nii.gz}_ants,${rootdir}/${vol_dir}/${vol//.nii.gz}_antsWarped.nii.gz,${rootdir}/${vol_dir}/${vol//.nii.gz}_antsInverseWarped.nii.gz] \
                --interpolation Linear \
                --winsorize-image-intensities [0.005,0.995] \
                --use-histogram-matching 1 \
                --transform Affine[0.1] \
                --metric CC[$T2_template,${vol%.nii.gz}_flirt.nii.gz,1,4] \
                --convergence [1000x350x250x100,1e-10,15] --shrink-factors 4x2x1x1 --smoothing-sigmas 3x2x1x0vox \
                --transform SyN[0.2,15,0] \
                --metric CC[$T2_template,${vol%.nii.gz}_flirt.nii.gz,1,4] \
                --convergence [150x100x75x50,1e-7,20] \
                --shrink-factors 4x2x1x1 \
                --smoothing-sigmas 1x0.5x0.25x0vox
            
            aff=${rootdir}/${vol_dir}/${vol//.nii.gz}_ants0GenericAffine.mat
            warp=${rootdir}/${vol_dir}/${vol//.nii.gz}_ants1Warp.nii.gz
            
        fi
        
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Applying registration"
        echo "  Using transforms:"
        echo "  - ${aff}"
        echo "  - ${warp}"
        flirt -in $vol -ref $T2_template -init ${vol_dir}_reg_output/iso_to_template.mat -applyxfm -out ${vol%.nii.gz}_reg_flirt.nii.gz
        antsApplyTransforms -d 3 \
            -i ${vol%.nii.gz}_reg_flirt.nii.gz \
            -r $T2_template \
            -t $aff \
            -t $warp \
            -o ${vol%.nii.gz}_reg_flirt_ants.nii.gz
            
        fslmaths ${vol%.nii.gz}_reg_flirt_ants.nii.gz -mas $post_reg_mask ${vol%.nii.gz}_reg_flirt_ants_masked.nii.gz
    #    fslmaths ${vol%.nii.gz}_reg_flirt.nii.gz -mas $post_reg_mask ${vol%.nii.gz}_reg_flirt_ants_masked.nii.gz
        
    done
    
    mv ${rootdir}/${vol_dir}/${vol_dir}_reg_output ${rootdir}/batch_reg_output/${vol_dir}_reg_output
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Making mirror halves"
    singularity exec $sing_path python /home/groups/rairan/PH_StrokeSC_MRI/NEW_AND_IMPROVED/sne2b-mirror.py ${rootdir}/batch_reg_output/ ${vol_dir}_reg_output
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] All done."
    
fi

mv ${vol_dir}_reg_output ../batch_reg_output
cd ..







