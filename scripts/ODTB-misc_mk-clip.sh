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

atlas_dir=/scratch/users/snewbank/atlases/

cd ${rootdir}/batch_output/batch_reg_output/${d}

mk=$(ls *_MK_reg.nii.gz)

if [ -f "${mk//.nii.gz}_orig.nii.gz" ]
then 
mv ${mk//.nii.gz}_orig.nii.gz $mk
rm ${mk//.nii.gz}_orig.nii.gz
fi

mv $mk ${mk//.nii.gz}_orig.nii.gz
fslmaths ${mk//.nii.gz}_orig.nii.gz -mas ${atlas_dir}/mouse/P56_MASK_100um.nii.gz ${mk//.nii.gz}_orig_masked.nii.gz

fslmaths ${mk//.nii.gz}_orig.nii.gz -max -5 -min 5 $mk

#p_hi=$(fslstats ${mk//.nii.gz}_orig_masked.nii.gz -P 99.5 | awk '{print $1}')
#p_lo=$(fslstats ${mk//.nii.gz}_orig_masked.nii.gz -P 0.5 | awk '{print $1}')
#fslmaths ${mk//.nii.gz}_orig.nii.gz -max $p_lo -min $p_hi $mk

rm ${mk//.nii.gz}_orig_masked.nii.gz