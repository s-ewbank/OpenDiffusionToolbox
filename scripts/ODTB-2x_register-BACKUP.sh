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
atlas_dir=/scratch/users/snewbank/atlases/

cd $rootdir
cd $vol_dir

#Set some variables
S0_vol=S0.nii.gz
FA_vol=FA.nii.gz
MD_vol=MD.nii.gz
pre_reg_mask=reg/pre_reg_mask.nii.gz
if [ ! -d "${vol_dir}_reg_output" ]; then mkdir ${vol_dir}_reg_output; fi

if [[ "$organism" == "mouse" ]]; then

    template=${atlas_dir}/mouse/P56_Atlas_100um.nii.gz
    init_txfm=${atlas_dir}/mouse/mouse_S0_to_ref.mat
    
elif [[ "$organism" == "rat" ]]; then

    template=${atlas_dir}/rat/WHS_SD_rat_BRAIN_T2.nii.gz
    init_txfm=${atlas_dir}/rat/rat_S0_to_ref.mat
    
elif [[ "$organism" == "human" ]]; then

    template=${atlas_dir}/human/FSL_HCP1065_FA_1mm.nii.gz
    
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Detected organism as ${organism}; registering to:"
echo "  - ${template}"


subj_txfm=reg/subj_to_template.mat

#########################################
# Skip all reg steps - just apply warps
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

mv ${vol_dir}_reg_output ../batch_reg_output
cd ..
