config_path=""

while [[ "$#" -gt 0 ]]; do
  case $1 in
    --config) config_path="$2"; shift ;;
    *) echo "Unknown parameter: $1" >&2; exit 1 ;;
  esac
  shift
done

if [[ -z "$config_path" || -z "$d" ]]; then
  echo "Usage: $0 --config <path/to/config.txt> --dir <directory_to_analyze>"
fi

filedir="$(cd "$(dirname "$BASH_SOURCE")" && pwd)"
source ${filedir}/../READ_CONFIG.sh $config_path

cd ${rootdir}/batch_output

cd ${rootdir}/batch_output
subdirs=$(ls)


echo "[$(date '+%Y-%m-%d %H:%M:%S')] Oops! Going back and computing RD"

count=1

for d in $subdirs
    do
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Computing RD for subdir $count called $d"
        cd $d
        singularity exec $container_path fslmaths MD.nii.gz -mul 3 -sub AD.nii.gz -div 2 RD.nii.gz
        cd ..
        count=$(($count+1))
    done