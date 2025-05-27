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
out_dir="${rootdir}/batch_tract_output/${out_fname}_tractography"

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

    label_path=${atlas_dir}/rat/PH_sub_atlas.nii.gz
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

if [ ! -d "${rootdir}/batch_tract_output" ]; then mkdir "${rootdir}/batch_tract_output"; fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] DSI script. Variables set as:"
echo "  - dir: ${d}"
echo "  - img: ${img}"
echo "  - aff: ${aff}"
echo "  - warp: ${warp}"


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

if ! [ -f ${out_dir}/${out_fname}_DSI.tt.gz.fullatlas_origspace.count.pass.connectogram.txt ]; then
    dsi_studio --action=ana \
        --source=${out_dir}/${out_fname}_DSI.fib.gz \
        --tract=${out_dir}/${out_fname}_DSI.tt.gz \
        --connectivity=${out_dir}/fullatlas_origspace.nii.gz \
        --thread_count=$SLURM_CPUS_PER_TASK \
        --connectivity_output=connectogram
fi



#########################################
# Doing special roi stats
#########################################

if [ "$special_ids" != "X" ]; then
    
    for special in "${special_ids[@]}"
    do
        if ! [ -f ${out_dir}/${out_fname}_${special}-filtered.tt.gz.stat.txt ]; then
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
            
            #dsi_studio --action=ana \
            #    --source=${out_dir}/${out_fname}_DSI.fib.gz \
            #    --tract=${out_dir}/${out_fname}_DSI.tt.gz \
            #    --roi=${out_dir}/atlas-roi-${special}_origspace.nii.gz \
            #    --thread_count=$SLURM_CPUS_PER_TASK \
            #    --output=${out_dir}/${out_fname}_${special}-filtered.tt.gz
            
            dsi_studio --action=ana \
                --source=${out_dir}/${out_fname}_DSI.fib.gz \
                --tract=${out_dir}/${out_fname}_${special}-filtered.tt.gz \
                --thread_count=$SLURM_CPUS_PER_TASK \
                --other_slices=${volume_paths_csv} \
                --export=stat
                
            #dsi_studio --action=ana \
            #    --source=${out_dir}/${out_fname}_DSI.fib.gz \
            #    --tract=${out_dir}/${out_fname}_${special}-filtered.tt.gz \
            #    --thread_count=$SLURM_CPUS_PER_TASK \
            #    --export=${out_dir}/${out_fname}_${special}-filtered.nii.gz
                
        fi
        
    done
    
    #if ! [ -f ${out_dir}/tract-stats.csv ]; then
    /opt/conda/envs/diffusionmritoolkit/bin/python ${scripts_dir}/ODTB-5b_dsi.py ${out_dir}
    #fi
fi