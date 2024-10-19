#!/bin/bash

#########################################
# 0 - GET ROOTDIR AND SCRIPTS DIR FROM CONFIG
#########################################
config_dir="$(dirname "$(realpath "$0")")"

rootdir=$(grep '^rootdir=' ${config_dir}/config.txt | cut -d'=' -f2)
scripts_dir=$(grep '^scripts_dir=' ${config_dir}/config.txt | cut -d'=' -f2)

mdt_container_path=$(grep '^mdt_container_path=' ${config_dir}/config.txt | cut -d'=' -f2)
miracl_container_path=$(grep '^miracl_container_path=' ${config_dir}/config.txt | cut -d'=' -f2)
dsi_container_path=$(grep '^dsi_container_path=' ${config_dir}/config.txt | cut -d'=' -f2)

mailtype=$(grep '^mailtype=' ${config_dir}/config.txt | cut -d'=' -f2)
mailuser=$(grep '^mailuser=' ${config_dir}/config.txt | cut -d'=' -f2)
slurm_out_path=$(grep '^slurm_out_path=' ${config_dir}/config.txt | cut -d'=' -f2)


if [ $# -eq 0 ]; then
    echo "No input provided for 'step'!"
    echo "Usage: $0 [fit|check_fit|combine_fit_outputs|reg|roi|combine_roi_outputs|vba_avg|vba_baseline|vba_delta|tract]"
    exit 1
fi

step=$1


for file in $(find "$scripts_dir" -type f); do
    sed -i "s|^#SBATCH --mail-type=.*|#SBATCH --mail-type=${mailtype}|" "$file"
    sed -i "s|^#SBATCH --mail-user=.*|#SBATCH --mail-user=${mailuser}|" "$file"
    sed -i "s|^#SBATCH --output=.*|#SBATCH --output=${slurm_out_path}|" "$file"
done

if [[ "$step" == "fit" ]]; then
    #########################################
    # 1 - DTI AND NODDI FITTING
    #########################################
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Doing DTI and NODDI fitting"
    cd $rootdir
    subdirs=$(ls)
    for d in $subdirs
        do
            echo $d
            sbatch ${scripts_dir}/sne1a-batch_dti_mdt.sbatch $rootdir $d $scripts_dir $mdt_container_path $miracl_container_path
        done
        
elif [[ "$step" == "check_fit" ]]; then
    #########################################
    # 1b - GET MISSED DIRECTORIES FROM INITIAL FITTING
    #########################################
    # Because sometimes py-opencl throws an error (device not found)
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Checking all directories for DTI/NODDI output..."
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
    rm */filt.nii.gz
    rm */filt_mppca.nii.gz
    rm */protocol.prtcl
    
    cd ${rootdir}/batch2
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Rerunning on batch that didn't run before..."
    subdirs=$(ls)
    for d in $subdirs
        do
            echo $d
            sbatch ${scripts_dir}/sne1a-batch_dti_mdt.sbatch ${rootdir}/batch2 $d $scripts_dir $mdt_container_path $miracl_container_path
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
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Identified $n_subdirs subdirectories in root directory for registration in ${rootdir}/batch_output"
    echo " "
    count=1
    cp ${scripts_dir}/sne2-register_MIRACL.sbatch ${rootdir}/batch_output/sne2-register_MIRACL.sbatch
    mkdir batch_miracl_output
    
    for d in $subdirs
    do
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running registration script on subdir $count called $d."
        sbatch sne2-register_MIRACL.sbatch ${rootdir}/batch_output $d 1 1 $miracl_container_path
        count=$(($count+1))
    done

elif [[ "$step" == "roi" ]]; then
    #########################################
    # 3a1 - ROI-BASED ANALYSIS
    #########################################
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Doing ROI-based analysis"
    output_dir="${rootdir}/batch_output/batch_miracl_output/roi_analysis_output/"
    cd $rootdir
    subdirs=$(ls)
    n_subdirs=`find . -mindepth 1 -maxdepth 1 -type d | wc -l`
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Identified $n_subdirs subdirectories in root directory for roi analysis in $rootdir"
    echo " "
    
    mkdir $output_dir
    mkdir "${output_dir}/roi_masks/"
    
    labels_path="${rootdir}/batch_output/batch_miracl_output/roi_analysis_output/annotation_hemi_combined_50um_parent-level_2.nii.gz"
    singularity exec $miracl_container_path cp /code/atlases/ara/annotation/annotation_hemi_combined_50um_parent-level_2.nii.gz $labels_path
    
    lut_path="${rootdir}/batch_output/batch_miracl_output/roi_analysis_output/ara_mouse_structure_graph_hemi_combined.csv"
    singularity exec $miracl_container_path cp /code/atlases/ara/ara_mouse_structure_graph_hemi_combined.csv $lut_path
    
    count=1
    for d in $subdirs
    do
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running roi_counter script on subdir $count called $d."
        sbatch ${scripts_dir}/sne3b-roi_counter.sbatch $rootdir $d $lut_path $labels_path $output_dir
        count=$(($count+1))
    done

elif [[ "$step" == "combine_roi_outputs" ]]; then
    #########################################
    # 3a2 - ROI-BASED ANALYSIS - COMBINE OUTPUTS
    #########################################
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Combining outputs from roi-based analysis"
    output_dir="${rootdir}/batch_output/batch_miracl_output/roi_analysis_output/"
    cd $rootdir
    
    singularity exec $miracl_container_path python ${scripts_dir}/sne3e-roi_out_combiner.py $output_dir
    find "$output_dir" -type f ! -name '*combined*' -print0 | xargs -0 -I {} rm {}

elif [[ "$step" == "vba_avg" || "$step" == "vba_baseline" || "$step" == "vba_delta" ]]; then
    #########################################
    # 3b - VOXEL-BASED ANALYSIS
    #########################################
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Voxel-based analysis"
    
    ml biology
    ml fsl
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Making directories etc for VBA output"
    
    miracl_output_dir=${rootdir}/batch_output/batch_miracl_output/
    
    cd $miracl_output_dir
    dir=vba/avg_volumes/
    tmap_dir="${miracl_output_dir}/vba/tmaps/"
    baseline_tmap_dir="${miracl_output_dir}/vba/baseline_tmaps/"
    mkdir -p vba
    mkdir -p $dir
    mkdir -p $tmap_dir
    mkdir -p $baseline_tmap_dir
    mkdir -p "${tmap_dir}/individual_diffs"
    labels_path="${miracl_output_dir}/${dir}/annotation_hemi_combined_50um_parent-level_2.nii.gz"
    mask_path="${miracl_output_dir}/${dir}/atlas_mask.nii.gz"
    
    if [[ "$step" == "vba_avg" ]]
    then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Voxel-based analysis - making average volumes"
        singularity exec $miracl_container_path cp /code/atlases/ara/annotation/annotation_hemi_combined_50um_parent-level_2.nii.gz $labels_path
        fslmaths $labels_path -bin $mask_path
        
        volume_types=("FA" "MD" "ODI" "NDI")
        
        for vol in "${volume_types[@]}"
        do
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Making an average of ${vol} images"
            mkdir "${dir}/${vol}"
            find "$miracl_output_dir" -type f -name "*_output_${vol}_reg.nii.gz" -print0 | xargs -0 -I {} cp {} "${dir}/${vol}"
            singularity exec $miracl_container_path python ${scripts_dir}/sne4b-avg_vol.py $dir $vol $mask_path
            rm -r "${dir}/${vol}"
        done
    fi
    
    if [[ "$step" == "vba_delta" ]]
    then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Voxel-based analysis - doing delta t-maps"
        volume_types=("FA" "MD" "ODI" "NDI")
        for vol in "${volume_types[@]}"
        do
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Making t-maps for ${vol}"
            singularity exec $miracl_container_path python ${scripts_dir}/sne4c-tmaps.py $miracl_output_dir $mask_path $tmap_dir "${miracl_output_dir}/vba/avg_volumes/" 1 0 ${vol} 0
        done
    fi
    
    if [[ "$step" == "vba_baseline" ]]
    then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Voxel-based analysis - doing baseline t-maps"
        
        sbatch ${scripts_dir}/sne4d-tmaps_baseline.sbatch ${scripts_dir} ${miracl_container_path} $volume_types $miracl_output_dir $mask_path $baseline_tmap_dir "${miracl_output_dir}/vba/avg_volumes/" 1 0 1

    fi

elif [[ "$step" == "tract" ]]; then
    #########################################
    # 3c - TRACTOGRAPHY-BASED ANALYSIS
    #########################################
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Doing tractography-based analysis"


else
    echo "Choose a valid option! [fit | check_fit | combine_fit_outputs | reg | roi | combine_roi_outputs | vba_avg | vba_baseline | vba_delta | tract]"
    exit 1
fi

