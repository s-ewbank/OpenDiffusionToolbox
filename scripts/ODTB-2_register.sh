#########################################
# Initialize and get info from config
#########################################

export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$SLURM_CPUS_PER_TASK

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running registration script"

config_path=""
d=""

while [[ "$#" -gt 0 ]]; do
  case $1 in
    --config) config_path="$2"; shift ;;
    --dir) vol_dir="$2"; shift ;;
    *) echo "Unknown parameter: $1" >&2; exit 1 ;;
  esac
  shift
done

if [[ -z "$config_path" || -z "$d" ]]; then
  echo "Usage: $0 --config <path/to/config.txt> --dir <directory_to_analyze>"
fi

filedir="$(cd "$(dirname "$BASH_SOURCE")" && pwd)"
source ${filedir}/../READ_CONFIG.sh $config_path

rootdir=${rootdir}/batch_output/

cd $rootdir
cd $vol_dir
pwd
mkdir ${vol_dir}_reg_output

#Set some variables
S0_vol=S0.nii.gz
FA_vol=FA.nii.gz
MD_vol=MD.nii.gz
mkdir reg
pre_reg_mask=reg/pre_reg_mask.nii.gz

if [[ "$organism" == "mouse" ]]; then

    template=${atlas_path}/mouse/P56_Atlas_100um.nii.gz
    init_txfm=${atlas_path}/mouse/mouse_S0_to_ref.mat
    
elif [[ "$organism" == "rat" ]]; then

    template=${atlas_path}/rat/WHS_SD_rat_BRAIN_T2.nii.gz
    init_txfm=${atlas_path}/rat/rat_S0_to_ref.mat
    
elif [[ "$organism" == "human" ]]; then

    template=${atlas_path}/human/FSL_HCP1065_FA_1mm.nii.gz
    
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Detected organism as ${organism}; registering to:"
echo "  - ${template}"

#########################################
# Brain extraction
#########################################
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting brain extraction."

if [[ "$organism" == "mouse" ]]; then
    
    mkdir ${rootdir}/${vol_dir}/temp
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Bias correcting with FAST"
    # FAST for bias correcting 
    fast -l 1 -I 20 -B $S0_vol
    rm ${S0_vol//.nii.gz}_seg.nii.gz ${S0_vol//.nii.gz}_pve*.nii.gz ${S0_vol//.nii.gz}_mixeltype.nii.gz
    
    fslmaths ${S0_vol//.nii.gz}_restore.nii.gz -nan ${S0_vol//.nii.gz}_restore_nan.nii.gz
    realvox=$(fslstats ${S0_vol//.nii.gz}_restore_nan.nii.gz -s | awk '{print $1}')
    rm ${S0_vol//.nii.gz}_restore_nan.nii.gz
    
    if (( $(echo "$realvox > 0" | bc -l) )); then
    	echo "[$(date '+%Y-%m-%d %H:%M:%S')] The FAST bias-corrected S0 image has valid data - continuing."
    else
    	echo "[$(date '+%Y-%m-%d %H:%M:%S')] The FAST bias-corrected S0 image is all NaNs - adjusting for bias by alternative method."
    	rm ${S0_vol//.nii.gz}_restore.nii.gz
    	p50=$(fslstats $S0_vol -P 50 | awk '{print $1}')
    	fslmaths $S0_vol -s 0.5 -max $S0_vol -s 0.2 -div $p50 ${S0_vol//.nii.gz}_bias.nii.gz
    	fslmaths $S0_vol -div ${S0_vol//.nii.gz}_bias.nii.gz ${S0_vol//.nii.gz}_intermed.nii.gz
    	p10=$(fslstats ${S0_vol//.nii.gz}_intermed.nii.gz -P 10 | awk '{print $1}')
    	p85=$(fslstats ${S0_vol//.nii.gz}_intermed.nii.gz -P 85 | awk '{print $1}')
    	fslmaths ${S0_vol//.nii.gz}_intermed.nii.gz -max $p10 -min $p85 -mul 3 -add $S0_vol -div 4 ${S0_vol//.nii.gz}_restore.nii.gz
    	rm ${S0_vol//.nii.gz}_bias.nii.gz ${S0_vol//.nii.gz}_intermed.nii.gz
    
    fi
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Extracting brain"
    # Take square root and do initial thresholding at 40th percentile
    fslmaths ${S0_vol//.nii.gz}_restore.nii.gz -sqrt ${rootdir}/${vol_dir}/temp/sqrt.nii.gz
    fslmaths ${rootdir}/${vol_dir}/temp/sqrt.nii.gz -thrP 40 ${rootdir}/${vol_dir}/temp/thresh.nii.gz
    
    # Binarize, smooth, threshold and binarize to get a mask excluding peripheral voxels which have survived thresholding, then smooth to create a cloud and multiply by original image
    fslmaths ${rootdir}/${vol_dir}/temp/thresh.nii.gz -bin -s 1.5 -thrP 80 -bin -s 1.5 ${rootdir}/${vol_dir}/temp/cloud.nii.gz
    fslmaths ${rootdir}/${vol_dir}/temp/cloud.nii.gz -mul ${S0_vol//.nii.gz}_restore.nii.gz -thrP 20 ${rootdir}/${vol_dir}/temp/cloud_intermed.nii.gz
    fslmaths ${rootdir}/${vol_dir}/temp/cloud_intermed.nii.gz -bin -ero -ero -ero -s 0.5 -thr 0.1 -bin -dilD -dilD ${rootdir}/${vol_dir}/temp/cloud_mask.nii.gz
    fslmaths ${S0_vol//.nii.gz}_restore.nii.gz -mas ${rootdir}/${vol_dir}/temp/cloud_mask.nii.gz -thrP 30 -bin -ero -ero -dilD -dilD $pre_reg_mask
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Masking non-brain tissue and ventricles"
    fslmaths $MD_vol -mas $pre_reg_mask -thrP 95 -bin -mul -1 -add 1 reg/ventricles.nii.gz
    fslmaths ${S0_vol//.nii.gz}_restore.nii.gz -mas $pre_reg_mask -mas reg/ventricles.nii.gz reg/S0_masked.nii.gz
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Clip at 90th percentile and apply POW!!"
    #p90=$(fslstats reg/S0_masked.nii.gz -P 90 | bc)  
    p90=$(fslstats reg/S0_masked.nii.gz -P 90 | awk '{print $1}')
    fslmaths reg/S0_masked.nii.gz -min $p90 -pow 2 reg/S0_masked_min_pow.nii.gz
    
    in_vol=reg/S0_masked_min_pow.nii.gz
    
    rm -r ${rootdir}/${vol_dir}/temp

elif [[ "$organism" == "rat" ]]; then
        
    
    mkdir ${rootdir}/${vol_dir}/temp
    
	echo "[$(date '+%Y-%m-%d %H:%M:%S')] Doing bias correction-like step for rat."
	p50=$(fslstats $S0_vol -P 50 | awk '{print $1}')
	fslmaths $S0_vol -s 0.5 -max $S0_vol -s 0.2 -div $p50 ${S0_vol//.nii.gz}_bias.nii.gz
	fslmaths $S0_vol -div ${S0_vol//.nii.gz}_bias.nii.gz ${S0_vol//.nii.gz}_intermed.nii.gz
	p10=$(fslstats ${S0_vol//.nii.gz}_intermed.nii.gz -P 10 | awk '{print $1}')
	p85=$(fslstats ${S0_vol//.nii.gz}_intermed.nii.gz -P 85 | awk '{print $1}')
	fslmaths ${S0_vol//.nii.gz}_intermed.nii.gz -max $p10 -min $p85 -mul 3 -add $S0_vol -div 4 ${S0_vol//.nii.gz}_restore.nii.gz
	rm ${S0_vol//.nii.gz}_bias.nii.gz ${S0_vol//.nii.gz}_intermed.nii.gz
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Extracting brain"
    # Take square root and do initial thresholding at 40th percentile
    fslmaths ${S0_vol//.nii.gz}_restore.nii.gz -sqrt ${rootdir}/${vol_dir}/temp/sqrt.nii.gz
    fslmaths ${rootdir}/${vol_dir}/temp/sqrt.nii.gz -thrP 40 ${rootdir}/${vol_dir}/temp/thresh.nii.gz
    
    p85=$(fslstats ${rootdir}/${vol_dir}/temp/thresh.nii.gz -P 85 | awk '{print $1}')
    fslmaths ${rootdir}/${vol_dir}/temp/thresh.nii.gz -min $p85 ${rootdir}/${vol_dir}/temp/thresh.nii.gz
    
    # Binarize, smooth, threshold and binarize to get a mask excluding peripheral voxels which have survived thresholding, then smooth to create a cloud and multiply by original image
    fslmaths ${rootdir}/${vol_dir}/temp/thresh.nii.gz -bin -s 4.5 -thrP 80 -bin -s 4.5 ${rootdir}/${vol_dir}/temp/cloud.nii.gz
    fslmaths ${rootdir}/${vol_dir}/temp/cloud.nii.gz -mul ${S0_vol//.nii.gz}_restore.nii.gz -thrP 20 ${rootdir}/${vol_dir}/temp/cloud_intermed.nii.gz
    fslmaths ${rootdir}/${vol_dir}/temp/cloud_intermed.nii.gz -bin -ero -ero -ero -s 1.5 -thr 0.1 -bin -dilD -dilD ${rootdir}/${vol_dir}/temp/cloud_mask.nii.gz
    fslmaths ${S0_vol//.nii.gz}_restore.nii.gz -mas ${rootdir}/${vol_dir}/temp/cloud_mask.nii.gz -thrP 20 -bin -ero -ero -dilD -dilD $pre_reg_mask
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Masking non-brain tissue"
    fslmaths ${S0_vol//.nii.gz}_restore.nii.gz -mas $pre_reg_mask reg/S0_masked.nii.gz
    
    fslmaths reg/S0_masked.nii.gz -bin -ero -ero -ero -s 1 \
        -mul reg/S0_masked.nii.gz -thrP 5 -bin -s 0.5 -thr 0.8 -bin -s 1 \
        -mul reg/S0_masked.nii.gz -thrP 5 -bin -s 0.5 -thr 0.95 -bin -s 1 \
        -mul reg/S0_masked.nii.gz -thrP 50 -bin -s 1.5 -thr 0.1 -bin ${rootdir}/${vol_dir}/temp/cloud_mask2.nii.gz
    
    fslmaths $pre_reg_mask -mas ${rootdir}/${vol_dir}/temp/cloud_mask2.nii.gz -s 0.5 -thr 0.5 -bin $pre_reg_mask
    fslmaths ${S0_vol//.nii.gz}_restore.nii.gz -mas $pre_reg_mask reg/S0_masked.nii.gz
    
    fslreorient2std $pre_reg_mask $pre_reg_mask
    fslreorient2std reg/S0_masked.nii.gz reg/S0_masked.nii.gz
    
    in_vol=reg/S0_masked.nii.gz
    
    rm -r ${rootdir}/${vol_dir}/temp


elif [[ "$organism" == "human" ]]; then

    fslmaths $S0_vol -bin $pre_reg_mask
    
    in_vol=$FA_vol
    
fi

#########################################
# FLIRT Registration
#########################################
echo "[$(date '+%Y-%m-%d %H:%M:%S')] FSL FLIRT registration step"

mv $in_vol reg

subj_txfm=reg/subj_to_template.mat

if [[ "$organism" == "rat" || "$organism" == "mouse" ]]; then

    flirt -in $in_vol \
        -ref $template \
        -init $init_txfm \
        -searchrx -15 15 -searchry -15 15 -searchrz -15 15 \
        -omat $subj_txfm \
        -noresampblur \
        -out ${in_vol//.nii.gz}_flirt.nii.gz
    
elif [[ "$organism" == "human" ]]; then

    flirt -in $in_vol \
        -ref $template \
        -omat $subj_txfm \
        -noresampblur \
        -out ${in_vol//.nii.gz}_flirt.nii.gz

fi

#########################################
# ANTs Registration
#########################################
echo "[$(date '+%Y-%m-%d %H:%M:%S')] ANTs registration step"

fslmaths $template -bin -dilD -dilD -dilD -dilD reg/template_mask.nii.gz

fslmaths ${in_vol//.nii.gz}_flirt.nii.gz -mas reg/template_mask.nii.gz ${in_vol//.nii.gz}_flirt.nii.gz

if [[ "$organism" == "mouse" ]]; then

    antsRegistration --dimensionality 3 --float 0 \
        --output [reg/subj_to_template,reg/subj_to_templateWarped.nii.gz,reg/subj_to_templateInverseWarped.nii.gz] \
        --interpolation Linear --use-histogram-matching 0 --winsorize-image-intensities [0.005,0.995] \
        --initial-moving-transform [${in_vol//.nii.gz}_flirt.nii.gz,$template,1] \
        --transform Affine[0.1] \
        --metric CC[${in_vol//.nii.gz}_flirt.nii.gz,$template,1,8] \
        --convergence [2000x700x550x200,1e-10,15] --shrink-factors 4x2x1x1 --smoothing-sigmas 3x2x1x0vox \
        --transform BSplineSyN[0.1,26,0,3] \
        --metric CC[${in_vol//.nii.gz}_flirt.nii.gz,$template,1,8] \
        --convergence [100x70x50x20,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox
        
elif [[ "$organism" == "rat" || "$organism" == "human" ]]; then

    antsRegistration --dimensionality 3 --float 0 \
        --output [reg/subj_to_template,reg/subj_to_templateWarped.nii.gz,reg/subj_to_templateInverseWarped.nii.gz] \
        --interpolation Linear --use-histogram-matching 0 --winsorize-image-intensities [0.005,0.995] \
        --initial-moving-transform [${in_vol//.nii.gz}_flirt.nii.gz,$template,1] \
        --transform Rigid[0.1] \
        --metric MI[${in_vol//.nii.gz}_flirt.nii.gz,$template,1,32,Regular,0.25] \
        --convergence [1000x500x250,1e-6,10] --shrink-factors 4x2x1 --smoothing-sigmas 3x2x1vox \
        --transform Affine[0.05] \
        --metric MI[${in_vol//.nii.gz}_flirt.nii.gz,$template,1,32,Regular,0.25] \
        --convergence [1000x500x250,1e-6,10] --shrink-factors 4x2x1 --smoothing-sigmas 3x2x1vox \
        --transform SyN[0.1,3,0] \
        --metric CC[${in_vol//.nii.gz}_flirt.nii.gz,$template,1,6] \
        --convergence [30x15x5,1e-6,10] --shrink-factors 4x2x1 --smoothing-sigmas 3x2x1vox

fi

#########################################
# Applying registration to other volumes
#########################################
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Applying registration to other volumes"

volumes=(ls *.nii.gz)

echo "Volumes to be registered are:"
for vol in "${volumes[@]}"
do
    echo "  - ${vol}"
done


for vol in "${volumes[@]}"
do
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Working on volume ${vol}"
    
    fslreorient2std $vol $vol
    
    flirt -in $vol \
        -ref $template \
        -init $subj_txfm \
        -applyxfm -noresampblur \
        -out ${vol//.nii.gz}_flirt.nii.gz
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Registering"
    antsApplyTransforms -d 3 \
        -r $template \
        -i ${vol//.nii.gz}_flirt.nii.gz -n Linear \
        -t [reg/subj_to_template0GenericAffine.mat,1] [reg/subj_to_template1InverseWarp.nii.gz,0] \
        -o ${vol%.nii.gz}_reg.nii.gz --float -v
    
    cp ${vol%.nii.gz}_reg.nii.gz ${vol_dir}_reg_output/${vol_dir}_${vol%.nii.gz}_reg.nii.gz
    
done

mv ${vol_dir}_reg_output ../batch_reg_output
cd ..

