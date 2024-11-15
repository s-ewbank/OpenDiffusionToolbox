#!/bin/bash

#########################################
# 0 - GET ROOTDIR AND SCRIPTS DIR FROM CONFIG
#########################################
config_path=/home/groups/rairan/NODDI/SCRIPTS2/config_mouse.txt

rootdir=$(grep '^rootdir=' ${config_path} | cut -d'=' -f2)
groups_file=$(grep '^groups_file=' ${config_path} | cut -d'=' -f2)

scripts_dir=$(grep '^scripts_dir=' ${config_path} | cut -d'=' -f2)

mdt_container_path=$(grep '^mdt_container_path=' ${config_path} | cut -d'=' -f2)
miracl_container_path=$(grep '^miracl_container_path=' ${config_path} | cut -d'=' -f2)
dsi_container_path=$(grep '^dsi_container_path=' ${config_path} | cut -d'=' -f2)

mailtype=$(grep '^mailtype=' ${config_path} | cut -d'=' -f2)
mailuser=$(grep '^mailuser=' ${config_path} | cut -d'=' -f2)
slurm_out_path=$(grep '^slurm_out_path=' ${config_path} | cut -d'=' -f2)

organism=$(grep '^organism=' ${config_path} | cut -d'=' -f2)
volumes=$(grep '^volumes=' ${config_path} | cut -d'=' -f2)
volumes=($(echo $volumes | tr ',' ' '))
do_mppca=$(grep '^do_mppca=' ${config_path} | cut -d'=' -f2)
do_slice=$(grep '^do_slice=' ${config_path} | cut -d'=' -f2)


if [ $# -eq 0 ]; then
    echo "No input provided for 'step'!"
    echo "Usage: $0 [fit|check_fit|combine_fit_outputs|reg|roi|combine_roi_outputs|vba_avg|vba_baseline|vba_delta|tract]"
    exit 1
fi

step=$1


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
    cd $rootdir
    mkdir batch_output
    subdirs=$(ls)
    for d in $subdirs
        do
            echo $d
            sbatch ${scripts_dir}/sne1a-batch_fit.sbatch $rootdir $d $scripts_dir $mdt_container_path $miracl_container_path $dsi_container_path $organism $do_mppca $do_slice "${volumes[@]}"
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
    rm */mppca.nii.gz
    rm */protocol.prtcl
    
    cd ${rootdir}/batch2
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Rerunning on batch that didn't run before..."
    subdirs=$(ls)
    for d in $subdirs
        do
            echo $d
            sbatch ${scripts_dir}/sne1a-batch_fit.sbatch ${rootdir}/batch2/ $d $scripts_dir $mdt_container_path $miracl_container_path $dsi_container_path $organism $do_mppca $do_slice "${volumes[@]}"
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
    cp ${scripts_dir}/sne2z-registerz.sbatch ${rootdir}/batch_output/sne2-register.sbatch
    mkdir batch_reg_output
    
    for d in $subdirs
    do
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running registration script on subdir $count called $d."
        sbatch sne2-register.sbatch $scripts_dir ${rootdir}/batch_output $d 1 1 $miracl_container_path $organism $volumes
        count=$(($count+1))
    done

elif [[ "$step" == "roi" ]]; then
    #########################################
    # 3a1 - ROI-BASED ANALYSIS
    #########################################
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Doing ROI-based analysis"
    ml biology
    ml fsl
    output_dir="${rootdir}/batch_output/batch_reg_output/roi_analysis_output/"
    cd ${rootdir}/batch_output/batch_reg_output/
    subdirs=$(ls)
    n_subdirs=`find . -mindepth 1 -maxdepth 1 -type d | wc -l`
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Identified $n_subdirs subdirectories in root directory for roi analysis in $rootdir"
    echo " "
    
    mkdir $output_dir
    
    labels_path="${rootdir}/batch_output/batch_reg_output/roi_analysis_output/labels.nii.gz"
    lut_path="${rootdir}/batch_output/batch_reg_output/roi_analysis_output/lut.csv"
    
    if [[ "$organism" == "mouse" ]]
    then
        singularity exec $miracl_container_path cp /code/atlases/ara/annotation/annotation_hemi_combined_50um_parent-level_2.nii.gz $labels_path
        singularity exec $miracl_container_path cp /code/atlases/ara/ara_mouse_structure_graph_hemi_combined.csv $lut_path
    fi
    if [[ "$organism" == "rat" ]]
    then
        singularity exec $miracl_container_path cp /code/atlases/ara/annotation/annotation_hemi_combined_50um_parent-level_2.nii.gz $labels_path
        singularity exec $miracl_container_path cp /code/atlases/ara/ara_mouse_structure_graph_hemi_combined.csv $lut_path
    fi
    if [[ "$organism" == "human" ]]
    then
        cp $FSLDIR/data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr25-2mm.nii.gz ${labels_path//.nii.gz}_cort.nii.gz
        cp $FSLDIR/data/atlases/HarvardOxford/HarvardOxford-sub-maxprob-thr25-2mm.nii.gz ${labels_path//.nii.gz}_subcort.nii.gz
        fslmaths ${labels_path//.nii.gz}_cort.nii.gz -bin -sub 1 -mul -1 ${labels_path//.nii.gz}_cort_invmask.nii.gz
        fslmaths ${labels_path//.nii.gz}_subcort.nii.gz -bin ${labels_path//.nii.gz}_subcort_mask.nii.gz
        fslmaths ${labels_path//.nii.gz}_subcort.nii.gz -add 100 -mas ${labels_path//.nii.gz}_cort_invmask.nii.gz -add ${labels_path//.nii.gz}_cort.nii.gz -mas ${labels_path//.nii.gz}_subcort_mask.nii.gz $labels_path
        rm ${labels_path//.nii.gz}_cort_invmask.nii.gz ${labels_path//.nii.gz}_subcort_mask.nii.gz
        cp $FSLDIR/data/atlases/HarvardOxford-Cortical.xml ${lut_path//.csv}_cort.xml
        cp $FSLDIR/data/atlases/HarvardOxford-Subcortical.xml ${lut_path//.csv}_subcort.xml
    fi
    count=1
    for d in $subdirs
    do
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running roi_counter script on subdir $count called $d."
        bash ${scripts_dir}/sne3a-roi_counter.sbatch ${rootdir}/batch_output/batch_reg_output/ $d $scripts_dir $miracl_container_path $lut_path $labels_path $output_dir $organism
        count=$(($count+1))
    done

elif [[ "$step" == "combine_roi_outputs" ]]; then
    #########################################
    # 3a2 - ROI-BASED ANALYSIS - COMBINE OUTPUTS
    #########################################
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Combining outputs from roi-based analysis"
    output_dir="${rootdir}/batch_output/batch_reg_output/roi_analysis_output/"
    cd $rootdir
    
    singularity exec $miracl_container_path python ${scripts_dir}/sne3c-roi_out_combiner.py $output_dir
    find "$output_dir" -type f ! -name '*combined*' -print0 | xargs -0 -I {} rm {}

elif [[ "$step" == "vba_avg" || "$step" == "vba_baseline" || "$step" == "vba_delta" ]]; then
    #########################################
    # 3b - VOXEL-BASED ANALYSIS
    #########################################
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Voxel-based analysis"
    
    ml biology
    ml fsl
    
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
    labels_path="${reg_output_dir}/${dir}/atlas_labels.nii.gz"
    mask_path="${reg_output_dir}/${dir}/atlas_mask.nii.gz"
        
    if [[ "$step" == "vba_avg" ]]
    then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Voxel-based analysis - making average volumes"
        if [[ "$organism" == "mouse" ]]
        then
            singularity exec $miracl_container_path cp /code/atlases/ara/annotation/annotation_hemi_combined_50um_parent-level_2.nii.gz $labels_path
        fi
        if [[ "$organism" == "rat" ]]
        then
            singularity exec $miracl_container_path cp /code/atlases/ara/annotation/annotation_hemi_combined_50um_parent-level_2.nii.gz $labels_path
        fi
        if [[ "$organism" == "human" ]]
        then
            cp $FSLDIR/data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr25-2mm.nii.gz ${labels_path//.nii.gz}_cort.nii.gz
            cp $FSLDIR/data/atlases/HarvardOxford/HarvardOxford-sub-maxprob-thr25-2mm.nii.gz ${labels_path//.nii.gz}_subcort.nii.gz
            fslmaths ${labels_path//.nii.gz}_cort.nii.gz -add ${labels_path//.nii.gz}_subcort.nii.gz -bin $mask_path
        fi
        fslmaths $labels_path -bin $mask_path
            
        
        for vol in "${volumes[@]}"
        do
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Making an average of ${vol} images"
            mkdir "${dir}/${vol}"
            find "$reg_output_dir" -type f -name "*_output_${vol}_reg.nii.gz" -print0 | xargs -0 -I {} cp {} "${dir}/${vol}"
            singularity exec $miracl_container_path python ${scripts_dir}/sne4b-avg_vol.py $dir $vol $mask_path $organism
            rm -r "${dir}/${vol}"
        done
    fi
    
    if [[ "$step" == "vba_delta" ]]
    then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Voxel-based analysis - doing delta t-maps"
        for vol in "${volumes[@]}"
        do
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Making t-maps for ${vol}"
            singularity exec $miracl_container_path python ${scripts_dir}/sne4c-tmaps.py $reg_output_dir $mask_path $tmap_dir "${reg_output_dir}/vba/avg_volumes/" $groups_file $organism 1 0 ${vol} 0
        done
    fi
    
    if [[ "$step" == "vba_baseline" ]]
    then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Voxel-based analysis - doing baseline t-maps"
        bash ${scripts_dir}/sne4d-tmaps_baseline.sbatch ${scripts_dir} ${miracl_container_path} $reg_output_dir $mask_path $baseline_tmap_dir "${reg_output_dir}/vba/avg_volumes/" $groups_file $organism 1 0 1 "${volumes[@]}"

    fi

elif [[ "$step" == "tract" ]]; then
    #########################################
    # 3c - TRACTOGRAPHY-BASED ANALYSIS
    #########################################
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Doing tractography-based analysis"

    bash ${scripts_dir}/sne5z-dsi.sbatch $miracl_container_path $dsi_container_path $rootdir ${rootdir}/batch_output $organism

else
    echo "Choose a valid option! [fit | check_fit | combine_fit_outputs | reg | roi | combine_roi_outputs | vba_avg | vba_baseline | vba_delta | tract]"
    exit 1
fi

