#!/bin/bash
#SBATCH -p normal
#SBATCH -c 1
#SBATCH --cpus-per-task=4
#SBATCH --time 12:00:00
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
        
        FA_template=/scratch/groups/rairan/NODDI/24-11-07/batch_output/batch_reg_output_V1/vba/avg_volumes/avg_FA_masked.nii.gz
        
        cp $FA_template ${rootdir}/${vol_dir}/FA_template.nii.gz
        fslmaths ${rootdir}/${vol_dir}/FA_template.nii.gz -subsamp2 ${rootdir}/${vol_dir}/FA_template_ss2.nii.gz
        
        fslmaths $FA_vol -mas ${rootdir}/${vol_dir}/pre_reg_mask.nii.gz ${FA_vol//.nii.gz}_bet.nii.gz
        
        flirt -in ${FA_vol//.nii.gz}_bet.nii.gz \
            -ref ${rootdir}/${vol_dir}/FA_template_ss2.nii.gz \
            -init ${script_dir_name}/mouse_S0_to_ref.mat \
            -searchrx -15 15 -searchry -15 15 -searchrz -15 15 \
            -omat ${rootdir}/${vol_dir}/subj_FA_to_ref.mat \
            -noresampblur \
            -out ${FA_vol//.nii.gz}_flirt.nii.gz

        antsRegistration --dimensionality 3 --float 0 \
            --output [${rootdir}/${vol_dir}/mr_allen_reg/FA_allen_mr_ants,${rootdir}/${vol_dir}/mr_allen_reg/FA_allen_mr_antsWarped.nii.gz,${rootdir}/${vol_dir}/mr_allen_reg/FA_allen_mr_antsInverseWarped.nii.gz] \
            --interpolation Linear --use-histogram-matching 0 --winsorize-image-intensities [0.005,0.995] \
            --initial-moving-transform [${FA_vol//.nii.gz}_flirt.nii.gz,${rootdir}/${vol_dir}/FA_template_ss2.nii.gz,1] \
            --transform Rigid[0.1] \
            --metric MI[${FA_vol//.nii.gz}_flirt.nii.gz,${rootdir}/${vol_dir}/FA_template_ss2.nii.gz,1,2] \
            --convergence [2000x700x550x200,1e-10,15] --shrink-factors 4x2x1x1 --smoothing-sigmas 3x2x1x0vox --transform Affine[0.1] \
            --metric CC[${FA_vol//.nii.gz}_flirt.nii.gz,${rootdir}/${vol_dir}/FA_template_ss2.nii.gz,1,8] \
            --convergence [2000x700x550x200,1e-10,15] --shrink-factors 4x2x1x1 --smoothing-sigmas 3x2x1x0vox --transform BSplineSyN[0.1,26,0,3] \
            --metric CC[${FA_vol//.nii.gz}_flirt.nii.gz,${rootdir}/${vol_dir}/FA_template_ss2.nii.gz,1,8] \
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
            
            flirt -in $vol \
                -ref ${rootdir}/${vol_dir}/FA_template.nii.gz \
                -init ${rootdir}/${vol_dir}/subj_FA_to_ref.mat \
                -applyxfm -noresampblur \
                -out ${vol//.nii.gz}_flirt.nii.gz
            
            flirt -in ${rootdir}/${vol_dir}/mr_allen_reg/FA_allen_mr_ants1InverseWarp.nii.gz -ref ${rootdir}/${vol_dir}/FA_template.nii.gz -applyxfm -out ${rootdir}/${vol_dir}/mr_allen_reg/FA_allen_mr_ants1InverseWarp_templateres.nii.gz
            
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Registering"
            singularity exec $miracl_container_path antsApplyTransforms -d 3 \
                -r ${rootdir}/${vol_dir}/FA_template.nii.gz \
                -i ${vol//.nii.gz}_flirt.nii.gz -n Linear \
                -t [${rootdir}/${vol_dir}/mr_allen_reg/FA_allen_mr_ants0GenericAffine.mat,1] [${rootdir}/${vol_dir}/mr_allen_reg/FA_allen_mr_ants1InverseWarp_templateres.nii.gz,0] \
                -o ${vol%.nii.gz}_reg.nii.gz --float -v
                
            rm $preproc_1
            rm $preproc_2
            rm $preproc_3
            rm $preproc_4
            #mv $preproc_4_masked ${vol%.nii.gz}_ss.nii.gz
            
        done
    
    fi
    
fi


mv ${vol_dir}_reg_output ../batch_reg_output
cd ..






