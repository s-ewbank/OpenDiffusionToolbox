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

rootdir=${rootdir}/batch_output/batch_reg_output/
output_dir=${rootdir}/roi_analysis_output/

mask_path=${rootdir}/${d}/atlas_mask.nii.gz

if [[ "$organism" == "mouse" ]]; then

    template_path=${atlas_path}/mouse/atlas_levels/ATLAS_ALL_100um.nii.gz
    lut_path=${atlas_path}/mouse/ABA_mouse_lut.csv
    cp ${atlas_path}/mouse/P56_MASK_100um.nii.gz $mask_path
    
elif [[ "$organism" == "rat" ]]; then

    template_path=${atlas_path}/rat/WHS_SD_rat_BRAIN_ATLAS.nii.gz
    lut_path=${atlas_path}/rat/WHS_SD_rat_lut.csv
    cp ${atlas_path}/rat/WHS_SD_rat_BRAIN_MASK.nii.gz $mask_path
    
elif [[ "$organism" == "human" ]]; then

    template_path=${atlas_path}/human/HarvardOxford-ALL-maxprob-thr50-1mm.nii.gz
    lut_path=${atlas_path}/human/human_lut.csv
    cp ${atlas_path}/human/atlas_mask.nii.gz $mask_path
    
fi


cd $rootdir

/opt/conda/envs/diffusionmritoolkit/bin/python ${scripts_dir}/ODTB-3b_roi-counter.py $d $lut_path $template_path $output_dir "${volumes[@]}"
