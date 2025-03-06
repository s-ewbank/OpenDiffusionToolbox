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
mkdir ${vol_dir}_reg_output

#Set some variables
S0_vol=S0.nii.gz
FA_vol=FA.nii.gz
MD_vol=MD.nii.gz

mkdir reg2
cp reg/pre_reg_mask.nii.gz reg2/pre_reg_mask.nii.gz
pre_reg_mask=reg2/pre_reg_mask.nii.gz
fslreorient2std $pre_reg_mask $pre_reg_mask

template=${rootdir}/batch_reg_output/vba/avg_volumes/avg_FA_masked.nii.gz
FLIRT_init=reg/subj_to_template.mat
ANTs_init_aff=reg/subj_to_template0GenericAffine.mat
ANTs_init_aff=reg/subj_to_template1InverseWarp.nii.gz

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Detected organism as ${organism}; registering to:"
echo "  - ${template}"

#########################################
# FLIRT Registration
#########################################
echo "[$(date '+%Y-%m-%d %H:%M:%S')] FSL FLIRT registration step"

fslmaths $FA_vol -mas $pre_reg_mask reg2/FA_masked.nii.gz
in_vol=reg2/FA_masked.nii.gz

subj_txfm=reg2/subj_to_template.mat

flirt -in $in_vol \
    -ref $template \
    -init $FLIRT_init \
    -searchrx -15 15 -searchry -15 15 -searchrz -15 15 \
    -omat $subj_txfm \
    -noresampblur \
    -out ${in_vol//.nii.gz}_flirt.nii.gz

#########################################
# ANTs Registration
#########################################
echo "[$(date '+%Y-%m-%d %H:%M:%S')] ANTs registration step"

fslmaths $template -bin -dilD -dilD -dilD -dilD reg2/template_mask.nii.gz

fslmaths ${in_vol//.nii.gz}_flirt.nii.gz -mas reg2/template_mask.nii.gz ${in_vol//.nii.gz}_flirt.nii.gz

if [[ "$organism" == "mouse" ]]; then

    antsRegistration --dimensionality 3 --float 0 \
        --output [reg2/subj_to_template,reg2/subj_to_templateWarped.nii.gz,reg2/subj_to_templateInverseWarped.nii.gz] \
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
        --output [reg2/subj_to_template,reg2/subj_to_templateWarped.nii.gz,reg2/subj_to_templateInverseWarped.nii.gz] \
        --interpolation Linear --use-histogram-matching 0 --winsorize-image-intensities [0.005,0.995] \
        --initial-moving-transform [${in_vol//.nii.gz}_flirt.nii.gz,$template,1] \
        --transform Affine[0.1] \
        --metric CC[${in_vol//.nii.gz}_flirt.nii.gz,$template,1,8] \
        --convergence [1000x500x250,1e-10,15] --shrink-factors 4x2x1 --smoothing-sigmas 2x1x0vox \
        --transform BSplineSyN[0.1,26,0,3] \
        --metric CC[${in_vol//.nii.gz}_flirt.nii.gz,$template,1,8] \
        --convergence [50x25x10,1e-6,10] --shrink-factors 4x2x1 --smoothing-sigmas 2x1x0vox
        
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

rm ${vol_dir}_reg_output/*_flirt_reg.nii.gz
rm ${vol_dir}_reg_output/*_reg_reg.nii.gz
mv ${vol_dir}_reg_output ../batch_reg_output2
cd ..
