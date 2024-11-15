#!/bin/bash
#SBATCH -p normal
#SBATCH -c 1
#SBATCH --cpus-per-task=4
#SBATCH --time 8:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=snewbank@stanford.edu
#SBATCH --output=/home/groups/rairan/NODDI/slurm_outputs//slurm-%j.out
#

ml biology
ml fsl
ml ants

miracl_sing_path=${1}
dsi_sing_path=${2}
rootdir=${3}
preregdir=${4}
organism=${5}

cd $rootdir
subdirs=$(ls)
n_subdirs=`find . -mindepth 1 -maxdepth 1 -type d | wc -l`
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Identified $n_subdirs subdirectories in root directory for DSI fiber tracking in $rootdir"
echo " "
count=1

if [[ "$organism" == "mouse" ]]
then
echo "bad!"
fi

if [[ "$organism" == "human" ]]
then
    for d in $subdirs; do
    
        d="${rootdir}/${d}/"
        out_fname=$(basename "$d")
        out_dir="${d}/${out_fname}_DSI"
        
        mkdir $out_dir
        
        omat=$(find "${preregdir}/batch_reg_output/" -type f -path "*/${out_fname}_output_reg_output/s0_to_ref.mat")
        ODI=$(find "$preregdir" -type f -path "*/${out_fname}_output/ODI.nii.gz")
        NDI=$(find "$preregdir" -type f -path "*/${out_fname}_output/NDI.nii.gz")
        
        img=`find ${d}/*.nii.gz`
        bval=`find ${d}/*.bval`
        bvec=`find ${d}/*.bvec`
        
        if [ ! -d "${rootdir}/batch_DSI_output" ]; then mkdir "${rootdir}/batch_DSI_output"; fi
        
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] DSI script. Variables set as:"
        echo "  - dir: ${d}"
        echo "  - img: ${img}"
        echo "  - omat: ${omat}"
        
        
        
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Making .src.gz"
        singularity exec $dsi_sing_path dsi_studio --action=src --source=$img --bval=$bval --bvec=$bvec --output=${out_dir}/${out_fname}_DSI.src.gz
        
        
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Applying flirt transforms"
        convert_xfm -omat ${out_dir}/inverted_transform.mat -inverse $omat
        flirt -in /scratch/users/snewbank/MKET2/data/batch_output/batch_reg_output/labels.nii.gz -ref $ODI -init ${out_dir}/inverted_transform.mat -applyxfm -noresampblur -out ${out_dir}/atlas_warp_flirt_RPI.nii.gz
        
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Reorienting and adjusting header"
        singularity exec $miracl_sing_path python /home/groups/rairan/NODDI/SCRIPTS2/scripts/sne5y-dsi.py $ODI ${out_dir}/atlas_warp_flirt_RPI.nii.gz ${out_dir}
        
        fslmaths ${out_dir}/atlas_warp_flirt_RPI.nii.gz -bin ${out_dir}/mask.nii.gz
        mv ${out_dir}/atlas_warp_abbrev.nii.gz ${out_dir}/atlas_warp_abbrev_uncleaned.nii.gz
        fslmaths ${out_dir}/atlas_warp_abbrev_uncleaned.nii.gz -bin -s 1 -thr 0.9 -bin ${out_dir}/atlas_warp_abbrev_mask.nii.gz
        fslmaths ${out_dir}/atlas_warp_abbrev_uncleaned.nii.gz -mas ${out_dir}/atlas_warp_abbrev_mask.nii.gz ${out_dir}/atlas_warp_abbrev.nii.gz
        
        singularity exec $dsi_sing_path dsi_studio --action=rec \
            --source=${out_dir}/${out_fname}_DSI.src.gz \
            --mask=${out_dir}/mask.nii.gz \
            --output=${out_dir}/${out_fname}_DSI.fib.gz
        
        singularity exec $dsi_sing_path dsi_studio --action=trk \
            --source=${out_dir}/${out_fname}_DSI.fib.gz \
            --output=${out_dir}/${out_fname}_DSI.tt.gz \
            --min_length=20 --max_length=300 --otsu_threshold=0.4 --fiber_count=1000000
        
        singularity exec $dsi_sing_path dsi_studio --action=ana \
            --source=${out_dir}/${out_fname}_DSI.fib.gz \
            --tract=${out_dir}/${out_fname}_DSI.tt.gz \
            --connectivity=${out_dir}/atlas_warp_abbrev.nii.gz \
            --other_slices=${ODI},${NDI}
            --connectivity_output=connectogram
            
        range=$(fslstats ${out_dir}/atlas_warp_abbrev.nii.gz -R)
        max=$(echo "$range" | awk '{print $2}')
        max_int=$(printf "%.0f" "$max")
        
        #for ((i=1; i<=max_int; i++)); do
        #    fslmaths ${out_dir}/atlas_warp_abbrev.nii.gz -thr ${i} -uthr ${i} -bin -s 0.1 -thr 0.75 -bin "${out_dir}/roi_${i}.nii.gz"
        #    singularity exec $dsi_sing_path dsi_studio --action=ana \
        #        --source=${out_dir}/${out_fname}_DSI.fib.gz \
        #        --tract=${out_dir}/${out_fname}_DSI.tt.gz \
        #        --roi=${out_dir}/roi_${i}.nii.gz \
        #        --export=stat --output=${out_dir}/${out_fname}_roi_${i}
        #    rm ${out_dir}/${out_fname}_roi_${i}.tt.gz
        #    rm "${out_dir}/roi_${i}.nii.gz"
        #done
        
        #singularity exec $dsi_sing_path dsi_studio --action=ana \
        #    --source=${out_dir}/${out_fname}_DSI.fib.gz \
        #    --tract=${out_dir}/${out_fname}_DSI.tt.gz \
        #    --regions=${out_dir}/atlas_warp_abbrev.nii.gz \
        #    --export=stat
        
        rm ${out_dir}/${out_fname}_DSI.src.gz
        mv $out_dir ${rootdir}/batch_DSI_output
        
    done
fi