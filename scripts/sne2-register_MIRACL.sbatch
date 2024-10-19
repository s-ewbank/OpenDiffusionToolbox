#!/bin/bash
#SBATCH -p normal
#SBATCH -c 1
#SBATCH --cpus-per-task=4
#SBATCH --time 48:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=snewbank@stanford.edu
#SBATCH --output=/home/groups/rairan/NODDI/slurm_outputs/slurm-%j.out
#

rootdir=${1}
vol_dir=${2}
do_ants=${3}
apply_ants_to_all_vols=${4}
miracl_container_path=${5}

echo "[$(date '+%Y-%m-%d %H:%M:%S')] MIRACL REGISTRATION SCRIPT"
echo " "

#######################################
# SET VARIABLES AND LOAD MODULES
#######################################
miracl_container_path=/scratch/groups/rairan/NODDI/MIRACL/miracl_latest.sif
labels=/code/atlases/ara/annotation/annotation_hemi_combined_50um.nii.gz
ml biology
ml fsl
 
#######################################
# MIRACL - SORTING OUT VOLUMES
#######################################
cd $rootdir
cd $vol_dir

mkdir mr_allen_reg
mkdir ${vol_dir}_miracl_output

cp FA.nii.gz ${vol_dir}_miracl_output/${vol_dir}_FA.nii.gz
cp MD.nii.gz ${vol_dir}_miracl_output/${vol_dir}_MD.nii.gz
cp NDI.nii.gz ${vol_dir}_miracl_output/${vol_dir}_NDI.nii.gz
cp ODI.nii.gz ${vol_dir}_miracl_output/${vol_dir}_ODI.nii.gz
t2_dat=S0.nii.gz

FA_vol=`find ${vol_dir}_miracl_output/*_FA.nii.gz`
MD_vol=`find ${vol_dir}_miracl_output/*_MD.nii.gz`
ODI_vol=`find ${vol_dir}_miracl_output/*_ODI.nii.gz`
NDI_vol=`find ${vol_dir}_miracl_output/*_NDI.nii.gz`
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

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running mri_allen_ants"
    singularity exec $miracl_container_path miracl reg mri_allen_ants -i ${t2_dat//.nii.gz}_grad.nii.gz \
        -o RSP -m combined -v 50 -b 0 -s 1 -f 0.3
        
fi

#######################################
# MIRACL - APPLYING ANTS TO OTHER VOLS
#######################################
if [[ $apply_ants_to_all_vols -ne 0 ]]
then
    mask=${rootdir}/${vol_dir}/${vol_dir}_miracl_output/miracl_mask.nii.gz
    preproc_1=${rootdir}/${vol_dir}/${vol_dir}_miracl_output/preproc1.nii.gz
    preproc_2=${rootdir}/${vol_dir}/${vol_dir}_miracl_output/preproc2.nii.gz
    preproc_3=${rootdir}/${vol_dir}/${vol_dir}_miracl_output/preproc3.nii.gz
    preproc_4=${rootdir}/${vol_dir}/${vol_dir}_miracl_output/preproc4.nii.gz
    preproc_4_masked=${rootdir}/${vol_dir}/${vol_dir}_miracl_output/preproc4_masked.nii.gz

    fslmaths ${rootdir}/${vol_dir}/mr_allen_reg/mr.nii.gz -bin $mask
    
    volumes=("$FA_vol" "$MD_vol" "$ODI_vol" "$NDI_vol")
    
    for vol in "${volumes[@]}"
    do
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Working on volume ${vol}"
        
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Doing reorient steps"
        singularity exec $miracl_container_path c3d $vol -orient RSP -o $preproc_1
        singularity exec $miracl_container_path c3d $preproc_1 -spacing 1.500000x1.495330x4.000000mm -o $preproc_2
        singularity exec $miracl_container_path c3d $preproc_2 -spacing .150000x.149533x.400000mm -o $preproc_3
        singularity exec $miracl_container_path flirt -in $preproc_3 -ref ${rootdir}/${vol_dir}/mr_allen_reg/mr.nii.gz -out $preproc_4 -applyxfm -usesqform 
        
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Masking"
        fslmaths $preproc_4 -mas $mask $preproc_4_masked
        
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Registering"
        singularity exec $miracl_container_path antsApplyTransforms -d 3 \
            -r /code/atlases/ara/template/average_template_50um_OBmasked.nii.gz \
            -i $preproc_4_masked -n Linear \
            -t [${rootdir}/${vol_dir}/mr_allen_reg/allen_mr_ants0GenericAffine.mat,1] [${rootdir}/${vol_dir}/mr_allen_reg/allen_mr_ants1InverseWarp.nii.gz,0] \
            -o ${vol%.nii.gz}_reg.nii.gz --float -v
            
        rm $preproc_1
        rm $preproc_2
        rm $preproc_3
        rm $preproc_4
        mv $preproc_4_masked ${vol%.nii.gz}_ss.nii.gz
        
    done

fi

mv ${vol_dir}_miracl_output ../batch_miracl_output
cd ..


        #echo "[$(date '+%Y-%m-%d %H:%M:%S')] Registering"
        #singularity exec $miracl_container_path antsApplyTransforms -d 3 \
        #    -r /code/atlases/ara/template/average_template_50um_OBmasked.nii.gz \
        #    -i $preproc_4_masked -n Linear \
        #    -t [${rootdir}/${vol_dir}/mr_allen_reg/allen_mr_ants1Warp.nii.gz,0] [${rootdir}/${vol_dir}/mr_allen_reg/allen_mr_ants0GenericAffine.mat,1] \
        #    -o out_w.nii.gz --float -v





