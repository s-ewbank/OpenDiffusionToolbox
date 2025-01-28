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

if [[ -z "$config_path" ]]; then
  echo "Usage: $0 --config <path/to/config.txt>"
fi

filedir="$(cd "$(dirname "$BASH_SOURCE")" && pwd)"
source ${filedir}/../READ_CONFIG.sh $config_path

#########################################
# TBSS
#########################################
regdir=${rootdir}/batch_output/batch_reg_output
tbss_dir=${regdir}/tbss
template=${regdir}/vba/avg_volumes/avg_FA_masked.nii.gz

if [ ! -d $tbss_dir ]; then mkdir $tbss_dir; fi

cd $d

mkdir stats

cp $template stats/mean_FA.nii.gz
fslmaths stats/mean_FA.nii.gz -bin stats/mean_FA_mask.nii.gz

cp *_FA_reg.nii.gz stats/all_FA.nii.gz

tbss_skeleton -i stats/mean_FA.nii.gz -o stats/mean_FA_skeleton.nii.gz
tbss_4_prestats 0.22

mkdir stats/${d}_tbss

for vol in "${volumes[@]}"
do
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Working on volume ${vol}"
    
    vol_reg=$(find . -type f -name "*_${vol}_reg.nii.gz")
    cp ${vol_reg} stats
    cd stats
    fslmaths ${vol_reg} -mas mean_FA_mask.nii.gz ${vol_reg}
    thresh=`cat thresh.txt`
    tbss_skeleton -i mean_FA -p $thresh mean_FA_skeleton_mask_dst ${FSLDIR}/data/standard/LowerCingulum_1mm all_FA ${d}_tbss/${vol_reg//_reg.nii.gz}_sk.nii.gz -a ${vol_reg}
    cd ..
    
done

mv stats/${d}_tbss $tbss_dir
rm -r stats