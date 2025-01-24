#!/bin/bash

#########################################
# 0 - GET USER INPUT THEN FIND ROOTDIR AND SCRIPTS DIR FROM CONFIG
#########################################

config_path=""
step=""

while [[ "$#" -gt 0 ]]; do
  case $1 in
    --config) config_path="$2"; shift ;;
    --step) step="$2"; shift ;;
    *) echo "Unknown parameter: $1" >&2; exit 1 ;;
  esac
  shift
done

if [[ -z "$config_path" || -z "$step" ]]; then
  echo "Usage: $0 --config <path/to/config.txt> --step <fit|check_fit|combine_fit_outputs|reg|roi|combine_roi_outputs|vba_avg|vba_baseline|vba_delta|tract> "
  exit 1
fi


filedir="$(cd "$(dirname "$BASH_SOURCE")" && pwd)"
source ${filedir}/READ_CONFIG.sh $config_path

for file in $(find "$scripts_dir" -type f); do
    sed -i "s|^#SBATCH --mail-type=.*|#SBATCH --mail-type=${mailtype}|" "$file"
    sed -i "s|^#SBATCH --mail-user=.*|#SBATCH --mail-user=${mailuser}|" "$file"
    sed -i "s|^#SBATCH --output=.*|#SBATCH --output=${slurm_out_path}/slurm-%j.out|" "$file"
done

if [[ "$step" == "fit" ]]; then
    #########################################
    # 1 - DTI AND NODDI FITTING
    #########################################
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Doing DTI and NODDI fitting"
    
    # Change job length
    sed -i 's/^#SBATCH --time .*/#SBATCH --time 6:00:00/' "${scripts_dir}/ODTB-0_job-submitter.sbatch"
    
    cd $rootdir
    mkdir batch_output
    subdirs=$(ls)
    
    for d in $subdirs
        do
            echo $d
            sbatch ${scripts_dir}/ODTB-0_job-submitter.sbatch --config $config_path --script ${scripts_dir}/ODTB-1a_fit.sh --dir $d
        done
        
elif [[ "$step" == "check_fit" ]]; then
    #########################################
    # 1b - GET MISSED DIRECTORIES FROM INITIAL FITTING
    #########################################
    # Because sometimes py-opencl throws an error (device not found)
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Checking all directories for DTI/NODDI output..."
    
    # Change job length
    sed -i 's/^#SBATCH --time .*/#SBATCH --time 6:00:00/' "${scripts_dir}/ODTB-0_job-submitter.sbatch"
    
    cd ${rootdir}/batch_output
    subdirs=$(ls)
    mkdir ../batch2
    
    for d in $subdirs
    do
    	if [ ! -f ${d}/ODI.nii.gz ]; then
    		echo "for ${d} - no NODDI output"
    		rm -r ${d}
    		cp -r ../${d//_output} ../batch2
    	fi
    done
    
    cd ../batch2
    rm -r */output
    rm */fullbinmask.nii.gz
    rm */raw_dwi_censored.nii.gz
    rm */raw_dwi_censored_mppca.nii.gz
    rm */tensorperfect_dwi.nii.gz
    
    cd ${rootdir}/batch2
    rootdir=${rootdir}/batch2
    subdirs=$(ls)
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Rerunning on batch that didn't run before..."
    
    for d in $subdirs
        do
            echo $d
            sbatch ${scripts_dir}/ODTB-0_job-submitter.sbatch --config $config_path --script ${scripts_dir}/ODTB-1a_fit.sh --dir $d
        done
        
elif [[ "$step" == "fit_backup" ]]; then
    #########################################
    # 1x - DTI AND NODDI FITTING - MDT only works here for now.
    #########################################
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Doing DTI and NODDI fitting"
    
    cd $rootdir
    mkdir batch_output
    subdirs=$(ls)
    n_subdirs=`find . -mindepth 1 -maxdepth 1 -type d | wc -l`
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Identified $n_subdirs subdirectories in root directory for initial fitting in ${rootdir}"
    count=1
    
    for d in $subdirs
        do
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running initial fitting script on subdir $count called $d."
            sbatch ${scripts_dir}/ODTB-1x_fit-BACKUP.sbatch --config $config_path --dir $d
            count=$(($count+1))
        done

elif [[ "$step" == "check_fit_backup" ]]; then
    #########################################
    # 1b - GET MISSED DIRECTORIES FROM INITIAL FITTING
    #########################################
    # Because sometimes py-opencl throws an error (device not found)
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Checking all directories for DTI/NODDI output..."
    
    # Change job length
    sed -i 's/^#SBATCH --time .*/#SBATCH --time 6:00:00/' "${scripts_dir}/ODTB-0_job-submitter.sbatch"
    
    cd ${rootdir}/batch_output
    subdirs=$(ls)
    mkdir ../batch2
    
    for d in $subdirs
    do
    	if [ ! -f ${d}/ODI.nii.gz ]; then
    		echo "for ${d} - no NODDI output"
    		rm -r ${d}
    		cp -r ../${d//_output} ../batch2
    	fi
    done
    
    cd ../batch2
    rm -r */output
    rm */fullbinmask.nii.gz
    rm */raw_dwi_censored.nii.gz
    rm */raw_dwi_censored_mppca.nii.gz
    rm */tensorperfect_dwi.nii.gz
    
    cd ${rootdir}/batch2
    rootdir=${rootdir}/batch2
    subdirs=$(ls)
    n_subdirs=`find . -mindepth 1 -maxdepth 1 -type d | wc -l`
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Identified $n_subdirs subdirectories from check which still need to be fitted in ${rootdir}"
    count=1
    
    for d in $subdirs
        do
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running initial fitting script on subdir $count called $d."
            sbatch ${scripts_dir}/ODTB-1x_fit-BACKUP.sbatch --config $config_path --dir $d --rerun 1
            count=$(($count+1))
        done

elif [[ "$step" == "combine_fit_outputs" ]]; then
    #########################################
    # 1c - COMBINE ORIGINAL FIT AND NEW FIT OUTPUTS
    #########################################
    # Because sometimes py-opencl throws an error (device not found)
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Combining outputs from fit and check_fit"
    
    mv ${rootdir}/batch2/batch_output/* ${rootdir}/batch_output/
    rm -r ${rootdir}/batch2/batch_output
    rm -r ${rootdir}/batch2

elif [[ "$step" == "reg" ]]; then
    #########################################
    # 2 - REGISTRATION
    #########################################
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Doing registration"
    
    cd ${rootdir}/batch_output
    subdirs=$(ls)
    n_subdirs=`find . -mindepth 1 -maxdepth 1 -type d | wc -l`
    
    # Change job length
    sed -i 's/^#SBATCH --time .*/#SBATCH --time 48:00:00/' "${scripts_dir}/ODTB-0_job-submitter.sbatch"
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Identified $n_subdirs subdirectories in root directory for registration in ${rootdir}/batch_output"
    echo " "
    count=1
    cp ${scripts_dir}/ODTB-2_register.sh ${rootdir}/batch_output/ODTB-2_register.sh
    mkdir batch_reg_output
    
    for d in $subdirs
    do
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running registration script on subdir $count called $d."
        sbatch ${scripts_dir}/ODTB-0_job-submitter.sbatch --config $config_path --script ${scripts_dir}/ODTB-2_register.sh --dir $d
        count=$(($count+1))
    done

elif [[ "$step" == "roi" ]]; then
    #########################################
    # 3a1 - ROI-BASED ANALYSIS
    #########################################
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Doing ROI-based analysis"
    
    postreg_root=${rootdir}/batch_output/batch_reg_output/
    output_dir="${postreg_root}/roi_analysis_output/"
    mkdir $output_dir
    cd ${postreg_root}
    subdirs=$(ls)
    n_subdirs=`find . -mindepth 1 -maxdepth 1 -type d | wc -l`
    
    # Change job length
    sed -i 's/^#SBATCH --time .*/#SBATCH --time 00:30:00/' "${scripts_dir}/ODTB-0_job-submitter.sbatch"
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Identified $n_subdirs subdirectories in root directory for roi analysis in $rootdir"
    echo " "
    
    count=1
    for d in $subdirs
    do
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running roi_counter script on subdir $count called $d."
        sbatch ${scripts_dir}/ODTB-0_job-submitter.sbatch --config $config_path --script ${scripts_dir}/ODTB-3a_roi-counter.sh --dir $d
        count=$(($count+1))
    done

elif [[ "$step" == "combine_roi_outputs" ]]; then
    #########################################
    # 3a2 - ROI-BASED ANALYSIS - COMBINE OUTPUTS
    #########################################
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Combining outputs from roi-based analysis"
    output_dir="${rootdir}/batch_output/batch_reg_output/roi_analysis_output/"
    cd $rootdir
    
    singularity exec $container_path /opt/conda/envs/diffusionmritoolkit/bin/python ${scripts_dir}/ODTB-3c_roi-out-combiner.py $output_dir "${volumes[@]}"
    find "$output_dir" -type f ! -name '*combined*' -print0 | xargs -0 -I {} rm {}

elif [[ "$step" == "vba_avg" || "$step" == "vba_baseline" || "$step" == "vba_delta" ]]; then
    #########################################
    # 3b - VOXEL-BASED ANALYSIS
    #########################################
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Voxel-based analysis"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Making directories etc for VBA output"
    
    reg_output_dir=${rootdir}/batch_output/batch_reg_output/
    
    cd $reg_output_dir
    dir=vba/avg_volumes/
    tmap_dir="${reg_output_dir}/vba/tmaps/"
    baseline_tmap_dir="${reg_output_dir}/vba/baseline_tmaps/"
    mkdir -p vba
    mkdir -p $dir
    mkdir -p $tmap_dir
    mkdir -p $baseline_tmap_dir
    mkdir -p "${tmap_dir}/individual_diffs"
    
    atlas_dir=/scratch/users/snewbank/atlases/
    labels_path="${reg_output_dir}/${dir}/atlas_labels.nii.gz"
    mask_path="${reg_output_dir}/${dir}/atlas_mask.nii.gz"
    
    if [[ "$step" == "vba_avg" ]]
    then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Voxel-based analysis - making average volumes"
        if [[ "$organism" == "mouse" ]]
        then
            cp ${atlas_dir}/mouse/P56_MASK_100um.nii.gz $mask_path
        fi
        if [[ "$organism" == "rat" ]]
        then
            cp ${atlas_dir}/rat/WHS_SD_rat_BRAIN_MASK.nii.gz $mask_path
        fi
        if [[ "$organism" == "human" ]]
        then
            fslmaths ${atlas_dir}/human/HarvardOxford-cort-maxprob-thr50-1mm.nii.gz \
                -add ${atlas_dir}/human/HarvardOxford-cort-maxprob-thr50-1mm.nii.gz \
                -bin $mask_path
        fi
        
        for vol in "${volumes[@]}"
        do
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Making an average of ${vol} images"
            mkdir "${dir}/${vol}"
            find "$reg_output_dir" -type f -name "*_output_${vol}_reg.nii.gz" -print0 | xargs -0 -I {} cp {} "${dir}/${vol}"
            singularity exec $container_path /opt/conda/envs/diffusionmritoolkit/bin/python ${scripts_dir}/ODTB-4b_avg-vol.py $dir $vol $mask_path $organism
            rm -r "${dir}/${vol}"
        done
    fi
    
    if [[ "$step" == "vba_delta" ]]
    then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Voxel-based analysis - doing delta t-maps"
        for vol in "${volumes[@]}"
        do
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Making t-maps for ${vol}"
            singularity exec $container_path /opt/conda/envs/diffusionmritoolkit/bin/python ${scripts_dir}/sne4c-tmaps.py $reg_output_dir $mask_path $tmap_dir "${reg_output_dir}/vba/avg_volumes/" $groups_file $organism 1 0 ${vol} 0
        done
    fi
    
    if [[ "$step" == "vba_baseline" ]]
    then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Voxel-based analysis - doing baseline t-maps"
        bash ${scripts_dir}/sne4d-tmaps_baseline.sbatch ${scripts_dir} ${miracl_container_path} $reg_output_dir $mask_path $baseline_tmap_dir "${reg_output_dir}/vba/avg_volumes/" $groups_file $organism 1 0 1 "${volumes[@]}"
    fi

elif [[ "$step" == "tract" ]]; then
    #########################################
    # 3c - TRACTOGRAPHY
    #########################################
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Doing tractography-based analysis"
    
    # Change job length
    sed -i 's/^#SBATCH --time .*/#SBATCH --time 0:30:00/' "${scripts_dir}/ODTB-0_job-submitter.sbatch"

    cd $rootdir
    subdirs=$(ls)
    n_subdirs=`find . -mindepth 1 -maxdepth 1 -type d | wc -l`
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Identified $n_subdirs subdirectories in root directory for DSI fiber tracking in $rootdir"
    echo " "
    count=1
    
    for d in $subdirs
    do
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running tractography script on subdir $count called $d."
        sbatch ${scripts_dir}/ODTB-0_job-submitter.sbatch --config $config_path --script ${scripts_dir}/ODTB-5a_dsi.sh --dir $d
        count=$(($count+1))
    done

elif [[ "$step" == "tbss" ]]; then
    #########################################
    # 3c - TRACT-BASED SPATIAL STATISTICS
    #########################################
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Doing tract-based spatial statistics"
    reg_output_dir=${rootdir}/batch_output/batch_reg_output/
    bash ${scripts_dir}/sne6a-tbss.sbatch $reg_output_dir

else
    echo "Choose a valid option! [fit | check_fit | combine_fit_outputs | reg | roi | combine_roi_outputs | vba_avg | vba_baseline | vba_delta | tract | tbss]"
    exit 1
fi

