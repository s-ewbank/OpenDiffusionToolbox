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
out_dir="${d}/${out_fname}_tractography"

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

if [ ! -d "${rootdir}/batch_tract_output" ]; then mkdir "${rootdir}/batch_tract_output"; fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] DSI script. Variables set as:"
echo "  - dir: ${d}"
echo "  - img: ${img}"
echo "  - aff: ${aff}"
echo "  - warp: ${warp}"


echo "[$(date '+%Y-%m-%d %H:%M:%S')] Making .src.gz"
dsi_studio --action=src --source=$img --bval=$bval --bvec=$bvec --output=${out_dir}/${out_fname}_DSI.src.gz

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Transforming label image and mask into original space"
convert_xfm -omat ${omat//.mat}_inverse.mat -inverse $omat
antsApplyTransforms -d 3 \
    -r $label_path \
    -i $label_path -n MultiLabel \
    -t [$aff,0] [$warp,0] \
    -o ${out_dir}/fullatlas_warp1.nii.gz --float -v
flirt -in ${out_dir}/fullatlas_warp1.nii.gz -ref ${FA} -init ${omat//.mat}_inverse.mat -applyxfm -usesqform -interp nearestneighbour -datatype short -out ${out_dir}/fullatlas_origspace.nii.gz
antsApplyTransforms -d 3 \
    -r $label_path \
    -i $mask -n MultiLabel \
    -t [$aff,0] [$warp,0] \
    -o ${out_dir}/mask_warp1.nii.gz --float -v
flirt -in ${out_dir}/mask_warp1.nii.gz -ref ${FA} -init ${omat//.mat}_inverse.mat -applyxfm -usesqform -interp nearestneighbour -out ${out_dir}/atlasmask_origspace.nii.gz
rm ${out_dir}/fullatlas_warp1.nii.gz ${out_dir}/mask_warp1.nii.gz
    
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Doing reconstruction and fiber tracking"
dsi_studio --action=rec \
    --source=${out_dir}/${out_fname}_DSI.src.gz \
    --mask=${out_dir}/atlasmask_origspace.nii.gz \
    --output=${out_dir}/${out_fname}_DSI.fib.gz

dsi_studio --action=trk \
    --source=${out_dir}/${out_fname}_DSI.fib.gz \
    --output=${out_dir}/${out_fname}_DSI.tt.gz \
    --min_length=2 --max_length=25 --tip_iteration=16 --fiber_count=2500000 --smoothing=0.5 --turning_angle=60 --method=1


#########################################
# Connectome
#########################################
volume_paths=()
for vol in "${volumes[@]}"; do
    volume_paths+=("${preregdir}/${out_fname}_output/${vol}.nii.gz")
done

IFS=,
volume_paths_csv="${volume_paths[*]}"
unset IFS

IFS=,
volume_names_csv="${volumes[*]}"
unset IFS

dsi_studio --action=ana \
    --source=${out_dir}/${out_fname}_DSI.fib.gz \
    --tract=${out_dir}/${out_fname}_DSI.tt.gz \
    --connectivity=${out_dir}/fullatlas_origspace.nii.gz \
    --other_slices=${volume_paths_csv} \
    --connectivity_value=${volume_names_csv} \
    --connectivity_output=connectogram


#########################################
# Doing special roi stats
#########################################
if [ "$special_ids" != "X" ]; then
    
    for special in "${special_ids[@]}"
    do
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Exporting stats for special roi with ID ${special}"
        if [ "$organism" != "rat" ]; then
            fslmaths $fulllabel_path -thr $special -uthr $special -Tmean -bin ${out_dir}/atlas-roi-${special}_atlasspace.nii.gz
            antsApplyTransforms -d 3 \
                -r $label_path \
                -i ${out_dir}/atlas-roi-${special}_atlasspace.nii.gz -n MultiLabel \
                -t [$aff,0] [$warp,0] \
                -o ${out_dir}/atlas-roi-${special}_warp1.nii.gz --float -v
            flirt -in ${out_dir}/atlas-roi-${special}_warp1.nii.gz -ref ${FA} -init ${omat//.mat}_inverse.mat -applyxfm -usesqform -interp nearestneighbour -out ${out_dir}/atlas-roi-${special}_origspace.nii.gz
            rm ${out_dir}/atlas-roi-${special}_atlasspace.nii.gz ${out_dir}/atlas-roi-${special}_warp1.nii.gz
        else
            fslmaths ${out_dir}/fullatlas_origspace.nii.gz -thr $special -uthr $special -bin ${out_dir}/atlas-roi-${special}_origspace.nii.gz
        fi
        
        dsi_studio --action=ana \
            --source=${out_dir}/${out_fname}_DSI.fib.gz \
            --tract=${out_dir}/${out_fname}_DSI.tt.gz \
            --roi=${out_dir}/atlas-roi-${special}_origspace.nii.gz \
            --output=${out_dir}/${out_fname}_${special}-filtered.tt.gz
        
        dsi_studio --action=ana \
            --source=${out_dir}/${out_fname}_DSI.fib.gz \
            --tract=${out_dir}/${out_fname}_${special}-filtered.tt.gz \
            --export=stat
            
        dsi_studio --action=ana \
            --source=${out_dir}/${out_fname}_DSI.fib.gz \
            --tract=${out_dir}/${out_fname}_${special}-filtered.tt.gz \
            --export=${out_dir}/${out_fname}_${special}-filtered.nii.gz
        
    done
    
    /opt/conda/envs/diffusionmritoolkit/bin/python ${scripts_dir}/ODTB-5b_dsi.py ${out_dir}
fi

mv $out_dir ${rootdir}/batch_tract_output
