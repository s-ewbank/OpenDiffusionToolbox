#!/bin/bash
#SBATCH -p normal
#SBATCH -c 1
#SBATCH --cpus-per-task=4
#SBATCH --time 6:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=snewbank@stanford.edu
#SBATCH --output=/scratch/users/snewbank/slurm_output//slurm-%j.out
#

config_path=""
d=""
rerun=0

while [[ "$#" -gt 0 ]]; do
  case $1 in
    --config) config_path="$2"; shift ;;
    --dir) d="$2"; shift ;;
    --rerun) rerun="$2"; shift ;;
    *) echo "Unknown parameter: $1" >&2; exit 1 ;;
  esac
  shift
done

if [[ -z "$config_path" || -z "$d" ]]; then
  echo "Usage: $0 --config <path/to/config.txt> --dir <directory_to_analyze>"
fi

filedir="$(cd "$(dirname "$BASH_SOURCE")" && pwd)"
source ${filedir}/../READ_CONFIG.sh $config_path

if [[ "$rerun" == 1 ]]; then
    $rootdir=${rootdir}/batch2
fi

dti=0
for volume in "${volumes[@]}"; do
    if [[ "$volume" == "FA" || "$volume" == "MD" || "$volume" == "AD" || "$volume" == "RD" ]]; then
        dti=1
        break
    fi
done

noddi=1
for volume in "${volumes[@]}"; do
    if [[ "$volume" == "NDI" || "$volume" == "ODI" ]]; then
        noddi=1
        break
    fi
done

dsi=0
for volume in "${volumes[@]}"; do
    if [[ "$volume" == "QA" ]]; then
        dsi=1
        break
    fi
done

dki=1
for volume in "${volumes[@]}"; do
    if [[ "$volume" == "MK" ]]; then
        dki=1
        break
    fi
done

#######################################
# SET VARIABLES
#######################################
protocol="protocol.prtcl"

cd $rootdir
ml system
ml clinfo

if [ ! -d "batch_output" ]; then mkdir batch_output; fi

#if [ ! -f "${rootdir}/filt_log.txt" ]; then touch "${rootdir}/filt_log.txt"; fi

subdirs=$(ls)
n_subdirs=`find . -mindepth 1 -maxdepth 1 -type d | wc -l`

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Variables set as:"
echo " -- Input data:"
echo "    > Root directory: $rootdir"
echo "    > Working directory: $d"
echo "    > Working directory: $d"
echo "    > Number subdirectories with presumed DWI volumes: $n_subdirs"
echo "    > List of subdirectories:"
for subdir in $subdirs; do echo "      + $subdir"; done
echo "    > Option for do_mppca: ${do_mppca}"
echo "    > Option for do_filt: ${do_slice}"
echo "    > Option for dki: ${dki}"

cd $d

if [ ! -d "output" ]; then mkdir output; fi
if [ ! -d "output/dti_output_dir" ]; then mkdir "output/dti_output_dir"; fi
if [ ! -d "output/noddi_output_dir" ] && [ "$noddi" -eq 1 ]; then mkdir "output/noddi_output_dir"; fi
if [ ! -d "output/kurtosis_output_dir" ] && [ "$dki" -eq 1 ]; then mkdir "output/kurtosis_output_dir"; fi

dwi=$(singularity exec $container_path find . -name "*.nii.gz" -exec sh -c 'for f; do [ "$(fslval "$f" dim4)" -gt 1 ] && echo "$f"; done' sh {} +)
bval=$(find . -name "*.bval")
bvec=$(find . -name "*.bvec")

#singularity exec $mdt_container_path env PYOPENCL_NO_CACHE=1 mdt-create-protocol $bvec $bval -o test.prtcl --cl-device-ind 0
#
#if [ ! -f "test.prtcl" ]; then
#    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Test protocol not created, so failure with pyopencl is inferred. Breaking script."
#    sbatch ${scripts_dir}/ODTB-1x_fit-with-mdt.sbatch --config $config_path --dir $d
#    exit 1
#fi
#
#echo "[$(date '+%Y-%m-%d %H:%M:%S')] Test protocol created, pyopencl must be working - proceeding with the script."
#rm test.prtcl

dipy_patch_rad=2
if [[ "$organism" == "mouse" ]]; then dipy_patch_rad=2; fi
if [[ "$organism" == "rat" ]]; then dipy_patch_rad=2; fi
if [[ "$organism" == "human" ]]; then dipy_patch_rad=2; fi

singularity exec $container_path fslmaths $dwi -Tmean -bin output/initial_dti_fullbinmask.nii.gz
singularity exec $container_path dtifit --data=$dwi --mask=output/initial_dti_fullbinmask.nii.gz --bvecs=$bvec --bvals=$bval --out=output/initial_dti

t2_dat="output/initial_dti_S0.nii.gz"
brainmask="output/initial_dti_BRAINMASK.nii.gz"
headmask="output/initial_dti_HEADMASK.nii.gz"

mkdir output/temp

#######################################
# MAKE BRAIN MASK
#######################################
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Making a brain mask"
if [[ "$organism" == "mouse" || "$organism" == "rat" ]]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Brain mask step 1 - Bias correcting"
    # FAST for bias correcting 
    
    singularity exec $container_path fslmaths $t2_dat -thrP 10 output/temp/masked.nii.gz
    p95=$(singularity exec $container_path fslstats output/temp/masked.nii.gz -P 95)
    singularity exec $container_path fslmaths output/temp/masked.nii.gz -min ${p95//.*} ${t2_dat//.nii.gz}_clip.nii.gz
    singularity exec $container_path fast -l 1 -I 20 -B -n 2 ${t2_dat//.nii.gz}_clip.nii.gz
    rm ${t2_dat//.nii.gz}_clip_seg.nii.gz ${t2_dat//.nii.gz}_clip_pve*.nii.gz ${t2_dat//.nii.gz}_clip_mixeltype.nii.gz
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Brain mask step 2 - Extracting brain"
    # Take square root and do initial thresholding at 40th percentile
    singularity exec $container_path fslmaths ${t2_dat//.nii.gz}_clip_restore.nii.gz -sqrt output/temp/sqrt.nii.gz
    singularity exec $container_path fslmaths output/temp/sqrt.nii.gz -thrP 40 output/temp/thresh.nii.gz
    
    # Binarize, smooth, threshold and binarize to get a mask excluding peripheral voxels which have survived thresholding, then smooth to create a cloud and multiply by original image
    singularity exec $container_path fslmaths output/temp/thresh.nii.gz -bin -s 1.5 -thrP 80 -bin -s 1.5 output/temp/cloud.nii.gz
    singularity exec $container_path fslmaths output/temp/cloud.nii.gz -mul ${t2_dat//.nii.gz}_clip_restore.nii.gz -thrP 20 output/temp/cloud_intermed.nii.gz
    singularity exec $container_path fslmaths output/temp/cloud_intermed.nii.gz -bin -ero -ero -ero -s 0.5 -thr 0.1 -bin -dilD -dilD output/temp/cloud_mask.nii.gz
    singularity exec $container_path fslmaths ${t2_dat//.nii.gz}_clip_restore.nii.gz -mas output/temp/cloud_mask.nii.gz -thrP 30 -bin -ero -ero -dilD -dilD $brainmask
    
    rm -r output/temp
else
    singularity exec $container_path bet $t2_dat ${t2_dat//.nii.gz}_bet.nii.gz -m
    mv ${t2_dat//.nii.gz}_bet_mask.nii.gz $brainmask
fi

singularity exec $container_path fslmaths $t2_dat -s 1 -thrP 15 -bin $headmask

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Initializing script for slice removal and MP-PCA"
singularity exec $container_path /opt/conda/envs/diffusionmritoolkit/bin/python ${scripts_dir}/ODTB-1b_filt-mppca.py ${rootdir}/${d} $dwi $bval $bvec output/initial_dti $do_mppca $do_slice $dipy_patch_rad

if [[ "$do_mppca" == 0 ]]; then
    if [[ "$do_slice" == 1 ]]; then
        dwi=$(find . -name "raw_dwi_censored.nii.gz")
    fi
else
    if [[ "$do_slice" == 0 ]]; then
        dwi=$(find . -name "raw_dwi_mppca.nii.gz")
    else
        dwi=$(find . -name "raw_dwi_censored_mppca.nii.gz")
    fi
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Switching to using ${dwi} as diffusion image after doing preprocessings steps."
    
#################################################################################
# EDDY CORRECT, MASK GENERATION, THEN DTIFIT
#################################################################################
if [[ "$do_slice" == 1 ]]; then
    for slice in slices/mask_*.nii.gz; do
        
        slice_i="${slice##*mask_}"
        slice_i="${slice_i%%.nii.gz*}"
        
        dwi=$(find slices/ -type f -name "dwi*_${slice_i}.nii.gz")
        mask=slices/mask_${slice_i}.nii.gz
        bval=slices/bval_${slice_i}.bval
        bvec=slices/bvec_${slice_i}.bvec
        protocol=slices/protocol_${slice_i}.prtcl
        
        #######################################
        # DTIFIT
        #######################################
        echo ""
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Fitting DTI to data with dtifit to slice ${slice_i}"
        echo ""
        singularity exec $container_path dtifit --data=$dwi --mask=$mask --bvecs=$bvec --bvals=$bval --out=output/dti_output_dir/dti
        
        #######################################
        # RUN MDT
        #######################################
        echo ""
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Fitting NODDI to data with MDT to slice ${slice_i}"
        echo ""
        singularity exec $mdt_container_path env PYOPENCL_NO_CACHE=1 mdt-create-protocol $bvec $bval -o $protocol --cl-device-ind 0
        singularity exec $mdt_container_path env PYOPENCL_NO_CACHE=1 mdt-model-fit -o output/noddi_output_dir/ "NODDI" $dwi $protocol $mask --cl-device-ind 0
        echo ""
        if [[ "$dki" == 1 ]]; then
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Fitting DKI to data with DIPY to slice ${slice_i}"
            echo ""
            singularity exec $container_path /opt/conda/envs/diffusionmritoolkit/bin/python ${scripts_dir}/ODTB-1c_dki.py output/kurtosis_output_dir/ $dwi $bval $bvec $mask
        fi
        echo ""
        
        if [ ! -d "output/final_output_dir" ]; then 
        mkdir "output/final_output_dir"
        cp output/dti_output_dir/dti_FA.nii.gz output/final_output_dir/FA.nii.gz
        cp output/dti_output_dir/dti_MD.nii.gz output/final_output_dir/MD.nii.gz
        cp output/dti_output_dir/dti_S0.nii.gz output/final_output_dir/S0.nii.gz
        cp output/dti_output_dir/dti_L1.nii.gz output/final_output_dir/AD.nii.gz
        singularity exec $container_path fslmaths output/dti_output_dir/dti_L2.nii.gz -add output/dti_output_dir/dti_L3.nii.gz -div 2 output/final_output_dir/RD.nii.gz
        cp output/noddi_output_dir/NODDI/ODI.nii.gz output/final_output_dir/ODI.nii.gz
        cp output/noddi_output_dir/NODDI/NDI.nii.gz output/final_output_dir/NDI.nii.gz
        if [[ "$dki" == 1 ]]; then
            cp output/kurtosis_output_dir/MK.nii.gz output/final_output_dir/MK.nii.gz
            rm -r output/kurtosis_output_dir
            mkdir output/kurtosis_output_dir
        fi
        rm -r output/noddi_output_dir
        rm -r output/dti_output_dir
        mkdir output/noddi_output_dir
        mkdir output/dti_output_dir
        else
        singularity exec $container_path fslmaths output/final_output_dir/FA.nii.gz -add output/dti_output_dir/dti_FA.nii.gz output/final_output_dir/FA.nii.gz
        singularity exec $container_path fslmaths output/final_output_dir/MD.nii.gz -add output/dti_output_dir/dti_MD.nii.gz output/final_output_dir/MD.nii.gz
        singularity exec $container_path fslmaths output/final_output_dir/S0.nii.gz -add output/dti_output_dir/dti_S0.nii.gz output/final_output_dir/S0.nii.gz
        singularity exec $container_path fslmaths output/final_output_dir/AD.nii.gz -add output/dti_output_dir/dti_L1.nii.gz output/final_output_dir/AD.nii.gz
        singularity exec $container_path fslmaths output/dti_output_dir/dti_L2.nii.gz -add output/dti_output_dir/dti_L3.nii.gz -div 2 -add output/final_output_dir/RD.nii.gz output/final_output_dir/RD.nii.gz
        singularity exec $container_path fslmaths output/final_output_dir/ODI.nii.gz -add output/noddi_output_dir/NODDI/ODI.nii.gz output/final_output_dir/ODI.nii.gz
        singularity exec $container_path fslmaths output/final_output_dir/NDI.nii.gz -add output/noddi_output_dir/NODDI/NDI.nii.gz output/final_output_dir/NDI.nii.gz
        if [[ "$dki" == 1 ]]; then
            singularity exec $container_path fslmaths output/final_output_dir/MK.nii.gz -add output/kurtosis_output_dir/MK.nii.gz output/final_output_dir/MK.nii.gz
            rm -r output/kurtosis_output_dir
            mkdir output/kurtosis_output_dir
        fi
        rm -r output/noddi_output_dir
        rm -r output/dti_output_dir
        mkdir output/noddi_output_dir
        mkdir output/dti_output_dir
        fi
        
    done
    
    mv slices/filt_log.csv output/final_output_dir/filt_log.csv
    rm -r slices
    
else
    
    mask=$headmask
    if [[ "$organism" == "mouse" || "$organism" == "rat" ]]; then
        mask=$headmask
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Detected organism as ${organism} and using mask ${headmask}"
    else
        mask=$brainmask
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Detected organism as ${organism} and using mask ${brainmask}"
    fi
    
    protocol=protocol.prtcl
    
    #######################################
    # DTIFIT
    #######################################
    echo ""
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Fitting DTI to data with dtifit"
    echo ""
    singularity exec $container_path dtifit --data=$dwi --mask=$mask --bvecs=$bvec --bvals=$bval --out=output/dti_output_dir/dti
    
    #######################################
    # RUN MDT
    #######################################
    if [[ "$noddi" == 1 ]]; then
        echo ""
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Fitting NODDI to data with MDT"
        echo ""
        singularity exec $mdt_container_path env PYOPENCL_NO_CACHE=1 mdt-create-protocol $bvec $bval -o $protocol
        singularity exec $mdt_container_path env PYOPENCL_NO_CACHE=1 mdt-model-fit -o output/noddi_output_dir/ "NODDI" $dwi $protocol $mask
    fi
    
    if [[ "$dki" == 1 ]]; then
        echo ""
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Fitting DKI to data with DIPY"
        echo ""
        singularity exec $container_path /opt/conda/envs/diffusionmritoolkit/bin/python ${scripts_dir}/ODTB-1c_dki.py output/kurtosis_output_dir/ $dwi $bval $bvec $mask
        echo ""
    fi
    
    mkdir "output/final_output_dir"
    cp output/dti_output_dir/dti_FA.nii.gz output/final_output_dir/FA.nii.gz
    cp output/dti_output_dir/dti_MD.nii.gz output/final_output_dir/MD.nii.gz
    cp output/dti_output_dir/dti_S0.nii.gz output/final_output_dir/S0.nii.gz
    cp output/dti_output_dir/dti_L1.nii.gz output/final_output_dir/AD.nii.gz
    singularity exec $container_path fslmaths output/dti_output_dir/dti_L2.nii.gz -add output/dti_output_dir/dti_L3.nii.gz -div 2 output/final_output_dir/RD.nii.gz
    rm -r output/dti_output_dir
    mkdir output/dti_output_dir
    if [[ "$noddi" == 1 ]]; then
        cp output/noddi_output_dir/NODDI/ODI.nii.gz output/final_output_dir/ODI.nii.gz
        cp output/noddi_output_dir/NODDI/NDI.nii.gz output/final_output_dir/NDI.nii.gz
        rm -r output/noddi_output_dir
        mkdir output/noddi_output_dir
    fi
    if [[ "$dki" == 1 ]]; then
        cp output/kurtosis_output_dir/MK.nii.gz output/final_output_dir/MK.nii.gz
        rm -r output/kurtosis_output_dir
        mkdir output/kurtosis_output_dir
    fi
        
fi


mv output/final_output_dir ${d}_output
mv ${d}_output ../batch_output/${d}_output






