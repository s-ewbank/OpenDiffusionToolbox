#########################################
# Initialize and get info from config
#########################################

export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$SLURM_CPUS_PER_TASK

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
out_dir=${rootdir}/batch_output/${out_fname}_output/origspace_template/
mkdir $out_dir

omat=$(find "$preregdir" -type f -path "*/${out_fname}_output/reg/subj_to_template.mat")
aff=$(find "$preregdir" -type f -path "*/${out_fname}_output/reg/subj_to_template0GenericAffine.mat")
warp=$(find "$preregdir" -type f -path "*/${out_fname}_output/reg/subj_to_template1InverseWarp.nii.gz")

FA=$(find "$preregdir" -type f -path "*/${out_fname}_output/FA.nii.gz")
bval=`find ${d}/*.bval`
bvec=`find ${d}/*.bvec`

if [ ! -f ${rootdir}/raw_means.txt ]; then
  touch ${rootdir}/raw_means.txt
fi

atlas_dir=/scratch/users/snewbank/atlases/

if [[ "$organism" == "mouse" ]]; then

    label_path=${atlas_dir}/mouse/atlas_levels/ATLAS_LVL6_100um.nii.gz
    fulllabel_path=${atlas_dir}/mouse/atlas_levels/ATLAS_ALL_100um.nii.gz
    mask=${rootdir}/mask.nii.gz
    #fslmaths ${atlas_dir}/mouse/P56_MASK_100um.nii.gz -ero -ero $mask
    
elif [[ "$organism" == "rat" ]]; then

    label_path=${atlas_dir}/rat/WHS_SD_rat_BRAIN_ATLAS.nii.gz
    mask=${atlas_dir}/rat/WHS_SD_rat_BRAIN_MASK.nii.gz
    
elif [[ "$organism" == "human" ]]; then

    template_path=${atlas_dir}/human/HarvardOxford-cort-maxprob-thr50-1mm.nii.gz
    fulllabel_path=${atlas_dir}/human/HarvardOxford-ALL-maxprob-thr50-1mm.nii.gz
    mask=${rootdir}/mask.nii.gz
    fslmaths ${atlas_dir}/human/HarvardOxford-cort-maxprob-thr50-1mm.nii.gz \
        -add ${atlas_dir}/human/HarvardOxford-cort-maxprob-thr50-1mm.nii.gz \
        -bin $mask
    
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Reverse reg volume quantification script. Variables set as:"
echo "  - dir: ${d}"
echo "  - aff: ${aff}"
echo "  - warp: ${warp}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Transforming label image and mask into original space"
convert_xfm -omat ${omat//.mat}_inverse.mat -inverse $omat
antsApplyTransforms -d 3 \
    -r $label_path \
    -i $mask -n MultiLabel \
    -t [$aff,0] [$warp,0] \
    -o ${out_dir}/mask_warp1.nii.gz --float -v
flirt -in ${out_dir}/mask_warp1.nii.gz -ref ${FA} -init ${omat//.mat}_inverse.mat -applyxfm -usesqform -interp nearestneighbour -out ${out_dir}/atlasmask_origspace.nii.gz
rm ${out_dir}/fullatlas_warp1.nii.gz ${out_dir}/mask_warp1.nii.gz

volume_paths=()
touch ${out_dir}/
for vol in "${volumes[@]}"; do
    fslmaths "${preregdir}/${out_fname}_output/${vol}.nii.gz" -mas ${out_dir}/atlasmask_origspace.nii.gz ${out_dir}/${vol}_masked.nii.gz
    stat=$(fslstats ${out_dir}/${vol}_masked.nii.gz -M)
    echo "$out_fname,$vol,$stat" >> ${rootdir}/raw_means.txt
done
