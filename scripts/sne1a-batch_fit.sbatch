#!/bin/bash
#SBATCH -p normal
#SBATCH -c 1
#SBATCH --cpus-per-task=4
#SBATCH --time 4:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=snewbank@stanford.edu
#SBATCH --output=/home/groups/rairan/NODDI/slurm_outputs//slurm-%j.out
#

rootdir=${1}
d=${2}
scripts_dir=${3}
mdt_container_path=${4}
miracl_container_path=${5}
dsi_container_path=${6}
organism=${7}
do_mppca=${8}
do_filt=${9}
volume_types=("${@:10}")

dti=0
for volume in "${volumes[@]}"; do
    if [[ "$volume" == "FA" || "$volume" == "MD" || "$volume" == "AD" || "$volume" == "RD" ]]; then
        dti=1
        break
    fi
done

noddi=0
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

#######################################
# SET VARIABLES
#######################################
protocol="protocol.prtcl"

cd $rootdir
ml biology
ml fsl
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
echo "    > Option for do_filt: ${do_filt}"

cd $d

if [ ! -d "output" ]; then mkdir output; fi
if [ ! -d "output/dti_output_dir" ]; then mkdir "output/dti_output_dir"; fi
if [ ! -d "output/noddi_output_dir" ]; then mkdir "output/noddi_output_dir"; fi

dwi=$(find . -name "*.nii.gz" -exec sh -c 'for f; do [ "$(fslval "$f" dim4)" -gt 1 ] && echo "$f"; done' sh {} +)
bval=$(find . -name "*.bval")
bvec=$(find . -name "*.bvec")

dipy_patch_rad=2
if [[ "$organism" == "mouse" ]]; then dipy_patch_rad=2; fi
if [[ "$organism" == "rat" ]]; then dipy_patch_rad=2; fi
if [[ "$organism" == "human" ]]; then dipy_patch_rad=16; fi

singularity exec $miracl_container_path python ${scripts_dir}/sne1b-filt-mppca.py ${rootdir}/${d} $dwi $bval $bvec $do_mppca $do_filt $dipy_patch_rad

if [[ "$do_mppca" == 0 ]]; then
    if [[ "$do_filt" == 1 ]]; then
        dwi=$(find . -name "filt.nii.gz")
    fi
else
    if [[ "$do_filt" == 0 ]]; then
        dwi=$(find . -name "mppca.nii.gz")
    else
        dwi=$(find . -name "filt_mppca.nii.gz")
    fi
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Switching to using ${dwi} as diffusion image after doing preprocessings steps."

if [[ "$organism" == "mouse" ]]; then
    
    #################################################################################
    # EDDY CORRECT, MASK GENERATION, THEN DTIFIT
    #################################################################################
    
    for slice in slices/filt_*.nii.gz; do
        
        slice_i="${slice##*filt_}"
        slice_i="${slice_i%%.nii.gz*}"
    
        mask=slices/mask_${slice_i}.nii.gz
        dwi=slices/filt_${slice_i}.nii.gz
        bval=slices/bval_${slice_i}.bval
        bvec=slices/bvec_${slice_i}.bvec
        protocol=slices/protocol_${slice_i}.prtcl
        
        #######################################
        # DTIFIT
        #######################################
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Fitting DTI to data with dtifit and outputting FA, MD, etc. in directory '${rootdir}/output/dti_output_dir/' with dtifit..."
        echo ""
        dtifit --data=$dwi --mask=$mask --bvecs=$bvec --bvals=$bval --kurt --kurtdir --out=output/dti_output_dir/dti
        
        #######################################
        # RUN MDT
        #######################################
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Done creating mask '$mask' and fitting DTI, now creating protocol $protocol in subdirectory '$d'..."
        echo ""
        singularity exec $mdt_container_path env PYOPENCL_NO_CACHE=1 mdt-create-protocol $bvec $bval -o $protocol
        
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Done creating '$protocol', now fitting NODDI to DWI volume $dwi..."
        singularity exec $mdt_container_path env PYOPENCL_NO_CACHE=1 mdt-model-fit -o output/noddi_output_dir/ "NODDI" $dwi $protocol $mask
        echo ""
        
        if [ ! -d "output/final_output_dir" ]; then 
        mkdir "output/final_output_dir"
        cp output/dti_output_dir/dti_FA.nii.gz output/final_output_dir/FA.nii.gz
        cp output/dti_output_dir/dti_MD.nii.gz output/final_output_dir/MD.nii.gz
        cp output/dti_output_dir/dti_S0.nii.gz output/final_output_dir/S0.nii.gz
        cp output/noddi_output_dir/NODDI/ODI.nii.gz output/final_output_dir/ODI.nii.gz
        cp output/noddi_output_dir/NODDI/NDI.nii.gz output/final_output_dir/NDI.nii.gz
        rm -r output/noddi_output_dir
        rm -r output/dti_output_dir
        mkdir output/noddi_output_dir
        mkdir output/dti_output_dir
        else
        fslmaths output/final_output_dir/FA.nii.gz -add output/dti_output_dir/dti_FA.nii.gz output/final_output_dir/FA.nii.gz
        fslmaths output/final_output_dir/MD.nii.gz -add output/dti_output_dir/dti_MD.nii.gz output/final_output_dir/MD.nii.gz
        fslmaths output/final_output_dir/S0.nii.gz -add output/dti_output_dir/dti_S0.nii.gz output/final_output_dir/S0.nii.gz
        fslmaths output/final_output_dir/ODI.nii.gz -add output/noddi_output_dir/NODDI/ODI.nii.gz output/final_output_dir/ODI.nii.gz
        fslmaths output/final_output_dir/NDI.nii.gz -add output/noddi_output_dir/NODDI/NDI.nii.gz output/final_output_dir/NDI.nii.gz
        rm -r output/noddi_output_dir
        rm -r output/dti_output_dir
        mkdir output/noddi_output_dir
        mkdir output/dti_output_dir
        fi
        
    done
    
    mv slices/filt_log.csv output/final_output_dir/filt_log.csv
    rm -r slices

fi

if [[ "$organism" == "human" ]]; then
        
    #################################################################################
    # EDDY CORRECT, MASK GENERATION, THEN DTIFIT
    #################################################################################
    mask=mask.nii.gz
    fslmaths $dwi -Tmean -s 5 -thrP 10 -bin $mask
        
    #######################################
    # DTIFIT
    #######################################
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Fitting DTI to data with dtifit and outputting FA, MD, etc. in directory '${rootdir}/output/dti_output_dir/' with dtifit..."
    echo ""
    dtifit --data=$dwi --mask=$mask --bvecs=$bvec --bvals=$bval --kurt --kurtdir --out=output/dti_output_dir/dti
    
    #######################################
    # RUN MDT
    #######################################
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Done creating mask '$mask' and fitting DTI, now creating protocol $protocol in subdirectory '$d'..."
    echo ""
    singularity exec $mdt_container_path env PYOPENCL_NO_CACHE=1 mdt-create-protocol $bvec $bval -o $protocol
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Done creating '$protocol', now fitting NODDI to DWI volume $dwi..."
    singularity exec $mdt_container_path env PYOPENCL_NO_CACHE=1 mdt-model-fit -o output/noddi_output_dir/ "NODDI" $dwi $protocol $mask
    echo ""
    
    mkdir "output/final_output_dir"
    cp output/dti_output_dir/dti_FA.nii.gz output/final_output_dir/FA.nii.gz
    cp output/dti_output_dir/dti_MD.nii.gz output/final_output_dir/MD.nii.gz
    cp output/dti_output_dir/dti_S0.nii.gz output/final_output_dir/S0.nii.gz
    cp output/noddi_output_dir/NODDI/ODI.nii.gz output/final_output_dir/ODI.nii.gz
    cp output/noddi_output_dir/NODDI/NDI.nii.gz output/final_output_dir/NDI.nii.gz
    rm -r output/noddi_output_dir
    rm -r output/dti_output_dir

fi

if [[ "$organism" == "rat" ]]; then
    
    threshP=25
    smooth=0.25
    kernel=1
    threshit=`fslstats $dwi -P $threshP`
    
    singularity exec $miracl_container_path python ${scripts_dir}/sne1b-filt-mppca.py ${rootdir}/${d}/ $dwi $bval $bvec 1 0
    dwi=${rootdir}/${d}/mppca.nii.gz
    
    fslmaths $dwi -s $smooth -thr $threshit -bin -kernel box ${kernel}x${kernel}x${kernel} -ero ${rootdir}/${d}/mask.nii.gz
    
    singularity exec $dsi_sing_path dsi_studio --action=src --source=$dwi --bval=$bval --bvec=$bvec --output=${dwi//.nii.gz}.src.gz
    singularity exec $dsi_sing_path dsi_studio --action=rec --source=${dwi//.nii.gz}.src.gz --mask=${rootdir}/${d}/mask.nii.gz --output=${dwi//.nii.gz}.fib.gz
    singularity exec $dsi_sing_path dsi_studio --action=ana --source=${dwi//.nii.gz}.fib.gz  --export=qa,iso,dti_fa,rd,ad,md
    
    mkdir ${rootdir}/${d}/${d}_output
    mv ${rootdir}/${d}/*.qa.nii.gz ${rootdir}/${d}/${d}_output/QA.nii.gz
    mv ${rootdir}/${d}/*.rd.nii.gz ${rootdir}/${d}/${d}_output/RD.nii.gz
    mv ${rootdir}/${d}/*.ad.nii.gz ${rootdir}/${d}/${d}_output/AD.nii.gz
    mv ${rootdir}/${d}/*.dti_fa.nii.gz ${rootdir}/${d}/${d}_output/FA.nii.gz
    mv ${rootdir}/${d}/*.iso.nii.gz ${rootdir}/${d}/${d}_output/ISO.nii.gz
    mv ${rootdir}/${d}/*.md.nii.gz ${rootdir}/${d}/${d}_output/MD.nii.gz
    
    mv ${d}/${subj}_output ${rootdir}/batch_output
    
fi


mv output/final_output_dir ${d}_output
mv ${d}_output ../batch_output/${d}_output







