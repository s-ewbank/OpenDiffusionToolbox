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

omat=$(find "$preregdir" -type f -path "*/${out_fname}_output/reg/subj_to_template.mat")
aff=$(find "$preregdir" -type f -path "*/${out_fname}_output/reg/subj_to_template0GenericAffine.mat")
warp=$(find "$preregdir" -type f -path "*/${out_fname}_output/reg/subj_to_template1InverseWarp.nii.gz")
FA=$(find "$preregdir" -type f -path "*/${out_fname}_output/FA.nii.gz")
reg_mask=$(find "$preregdir" -type f -path "*/masks/${out_fname}_output_FA_reg_mask.nii.gz")

if ! [ -f "${preregdir}/${out_fname}_output/postreg_mask_origspace.nii.gz" ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Transforming label image and mask into original space"
    convert_xfm -omat ${omat//.mat}_inverse.mat -inverse $omat
    antsApplyTransforms -d 3 \
        -r $reg_mask \
        -i $reg_mask -n MultiLabel \
        -t [$aff,0] [$warp,0] \
        -o ${preregdir}/${out_fname}_output/postreg_mask_warp1.nii.gz --float -v
    flirt -in ${preregdir}/${out_fname}_output/postreg_mask_warp1.nii.gz \
        -ref ${FA} -init ${omat//.mat}_inverse.mat -applyxfm -usesqform \
        -interp nearestneighbour -datatype short \
        -out ${preregdir}/${out_fname}_output/postreg_mask_origspace.nii.gz
    rm ${preregdir}/${out_fname}_output/postreg_mask_warp1.nii.gz
fi