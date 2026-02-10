#########################################
# Initialize and get info from config
#########################################

config_path=""
vol=""

while [[ "$#" -gt 0 ]]; do
  case $1 in
    --config) config_path="$2"; shift ;;
    --dir) vol="$2"; shift ;;
    *) echo "Unknown parameter: $1" >&2; exit 1 ;;
  esac
  shift
done

if [[ -z "$config_path" || -z "$d" ]]; then
  echo "Usage: $0 --config <path/to/config.txt> --dir <directory_to_analyze>"
fi

filedir="$(cd "$(dirname "$BASH_SOURCE")" && pwd)"
source ${filedir}/../READ_CONFIG.sh $config_path

#########################################
# Perform randomise
#########################################

reg_output_dir=${rootdir}/batch_output/batch_reg_output/

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Making design and contrast vest files"
save_stem=$(/opt/conda/envs/diffusionmritoolkit/bin/python ${scripts_dir}/ODTB-misc_make_design.py $groups_file ${reg_output_dir}/tbss/${vol}_)
Text2Vest ${reg_output_dir}/tbss/${vol}_${save_stem}_design.mat ${reg_output_dir}/tbss/${vol}_${save_stem}_design.mat
Text2Vest ${reg_output_dir}/tbss/${vol}_${save_stem}_design.con ${reg_output_dir}/tbss/${vol}_${save_stem}_design.con
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Made design and contrast vest files - now saving everything with comparison title ${vol}_${save_stem}"

randomise_inputs=()

cd ${reg_output_dir}
cd tbss

while IFS=',' read -r subject group timepoint filename; do
    filename=$(echo "$filename" | tr -d '\r')
    full_path=${filename}_output_reg_output_tbss/${filename}_output_${vol}_sk.nii.gz
    randomise_inputs+=("$full_path")
done < <(tail -n +2 "$groups_file")

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Making merged volume"
fslmerge -t ${reg_output_dir}/tbss/merged_${vol}_${save_stem}.nii.gz "${randomise_inputs[@]}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running randomise"
randomise -i ${reg_output_dir}/tbss/merged_${vol}_${save_stem}.nii.gz \
    -o ${reg_output_dir}/tbss/${vol}_${save_stem} \
    -d ${reg_output_dir}/tbss/${vol}_${save_stem}_design.mat \
    -t ${reg_output_dir}/tbss/${vol}_${save_stem}_design.con \
    -m ${reg_output_dir}/tbss/sk.nii.gz \
    -n 10000 --T2
    
mkdir ${vol}_${save_stem}
mv *${vol}_${save_stem}* ${vol}_${save_stem}
mkdir -p randomise_output
mv ${vol}_${save_stem} randomise_output




