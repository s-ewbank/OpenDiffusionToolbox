#########################################
# Initialize and get info from config
#########################################

config_path=""
d=""

while [[ "$#" -gt 0 ]]; do
  case $1 in
    --config) config_path="$2"; shift ;;
    --dir) d="$2"; shift ;;
    *) echo "Unknown parameter: $1" >&2; exit 1 ;;
  esac
  shift
done

if [[ -z "$config_path" || -z "$d" ]]; then
  echo "Usage: $0 --config <path/to/config.txt> --dir <directory_to_analyze>"
fi

filedir="$(cd "$(dirname "$BASH_SOURCE")" && pwd)"
source ${filedir}/../READ_CONFIG.sh $config_path

preregdir=${rootdir}/batch_output/

cd $rootdir
d="${rootdir}/${d}/"
out_fname=$(basename "$d")
out_dir="${d}/${out_fname}_prob_tractography"

mkdir $out_dir

omat=$(find "$preregdir" -type f -path "*/${out_fname}_output/reg/subj_to_template.mat")
aff=$(find "$preregdir" -type f -path "*/${out_fname}_output/reg/subj_to_template0GenericAffine.mat")
warp=$(find "$preregdir" -type f -path "*/${out_fname}_output/reg/subj_to_template1InverseWarp.nii.gz")

FA=$(find "$preregdir" -type f -path "*/${out_fname}_output/FA.nii.gz")
bval=`find ${d}/*.bval`
bvec=`find ${d}/*.bvec`

if [[ "$do_mppca" == 0 ]]; then
    if [[ "$do_slice" == 1 ]]; then
        img=$(find ${d} -name "raw_dwi_censored.nii.gz")
    fi
else
    if [[ "$do_slice" == 0 ]]; then
        img=$(find ${d} -name "raw_dwi_mppca.nii.gz")
    else
        img=$(find ${d} -name "raw_dwi_censored_mppca.nii.gz")
    fi
fi

atlas_dir=/scratch/users/snewbank/atlases/

if [[ "$organism" == "mouse" ]]; then

    label_path=${atlas_dir}/mouse/atlas_levels/ATLAS_LVL6_100um.nii.gz
    fulllabel_path=${atlas_dir}/mouse/atlas_levels/ATLAS_ALL_100um.nii.gz
    mask=${rootdir}/mask.nii.gz
    fslmaths ${atlas_dir}/mouse/P56_MASK_100um.nii.gz -ero -ero $mask
    
elif [[ "$organism" == "rat" ]]; then

    label_path=${atlas_dir}/rat/WHS_SD_rat_BRAIN_ATLAS.nii.gz
    fulllabel_path=${atlas_dir}/rat/WHS_SD_rat_BRAIN_ATLAS.nii.gz
    mask=${atlas_dir}/rat/WHS_SD_rat_BRAIN_MASK.nii.gz
    
elif [[ "$organism" == "human" ]]; then

    template_path=${atlas_dir}/human/HarvardOxford-cort-maxprob-thr50-1mm.nii.gz
    fulllabel_path=${atlas_dir}/human/HarvardOxford-ALL-maxprob-thr50-1mm.nii.gz
    mask=${rootdir}/mask.nii.gz
    fslmaths ${atlas_dir}/human/HarvardOxford-cort-maxprob-thr50-1mm.nii.gz \
        -add ${atlas_dir}/human/HarvardOxford-cort-maxprob-thr50-1mm.nii.gz \
        -bin $mask
    
fi

#########################################
# Fiber tracking
#########################################

if [ ! -d "${rootdir}/batch_prob_tract_output" ]; then mkdir "${rootdir}/batch_prob_tract_output"; fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] DSI script. Variables set as:"
echo "  - dir: ${d}"
echo "  - img: ${img}"
echo "  - aff: ${aff}"
echo "  - warp: ${warp}"



echo "[$(date '+%Y-%m-%d %H:%M:%S')] Transforming label image and mask into original space"
convert_xfm -omat ${omat//.mat}_inverse.mat -inverse $omat

antsApplyTransforms -d 3 \
    -r $mask \
    -i $mask -n MultiLabel \
    -t [$aff,0] [$warp,0] \
    -o ${out_dir}/mask_warp1.nii.gz --float -v
flirt -in ${out_dir}/mask_warp1.nii.gz -ref ${FA} -init ${omat//.mat}_inverse.mat -applyxfm -usesqform -interp nearestneighbour -out ${out_dir}/nodif_brain_mask.nii.gz

mean_1=$(fslstats ${out_dir}/mask_warp1.nii.gz -m)
mean_final=$(fslstats ${out_dir}/nodif_brain_mask.nii.gz -m)

if (( $(echo "$mean_1 <= 0" | bc -l) )); then
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] First mask warp is empty - exiting :("
  exit 1
fi

if (( $(echo "$mean_final <= 0" | bc -l) )); then
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Mask is empty - exiting :("
  exit 1
fi

rm ${out_dir}/mask_warp1.nii.gz

cp $img ${out_dir}/data.nii.gz
cp $bvec ${out_dir}/bvecs
cp $bval ${out_dir}/bvals

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running bedpostx"
bedpostx ${out_dir}

seed_vol_paths=""

for special in "${special_ids[@]}"
do
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Making ROI ID ${special}"
    fslmaths $fulllabel_path -thr $special -uthr $special -Tmean -bin ${out_dir}/atlas-roi-${special}_atlasspace.nii.gz
    antsApplyTransforms -d 3 \
        -r $mask \
        -i ${out_dir}/atlas-roi-${special}_atlasspace.nii.gz -n MultiLabel \
        -t [$aff,0] [$warp,0] \
        -o ${out_dir}/atlas-roi-${special}_warp1.nii.gz --float -v
    flirt -in ${out_dir}/atlas-roi-${special}_warp1.nii.gz -ref ${FA} -init ${omat//.mat}_inverse.mat -applyxfm -usesqform -interp nearestneighbour -out ${out_dir}/atlas-roi-${special}_origspace.nii.gz
    rm ${out_dir}/atlas-roi-${special}_atlasspace.nii.gz ${out_dir}/atlas-roi-${special}_warp1.nii.gz
    probtrackx2 -s ${out_dir}.bedpostX/merged -m ${out_dir}/nodif_brain_mask.nii.gz -x ${out_dir}/atlas-roi-${special}_origspace.nii.gz --dir=${out_dir}/probtrackx_${special} --ompl --opd
done

mkdir ${out_dir}/bedpostx
mv ${out_dir}.bedpostX/* ${out_dir}/bedpostx
rmdir ${out_dir}.bedpostX

mv $out_dir ${rootdir}/batch_prob_tract_output






