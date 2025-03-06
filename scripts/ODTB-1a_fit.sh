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

dki=0
for volume in "${volumes[@]}"; do
    if [[ "$volume" == "MK" ]]; then
        dki=1
        break
    fi
done

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Detected and will do the following:"
echo "  - do DTI: $dti"
echo "  - do NODDI: $noddi"
echo "  - do DKI: $dki"
echo "  - do DSI: $dsi"

#######################################
# SET VARIABLES
#######################################
protocol="protocol.prtcl"

cd $rootdir

if [ ! -d "batch_output" ]; then mkdir batch_output; fi

subdirs=$(ls)
n_subdirs=`find . -mindepth 1 -maxdepth 1 -type d | wc -l`

cd $d

if [ ! -d "output" ]; then mkdir output; fi
if [ ! -d "output/dti_output_dir" ]; then mkdir "output/dti_output_dir"; fi
if [ ! -d "output/noddi_output_dir" ] && [ "$noddi" -eq 1 ]; then mkdir "output/noddi_output_dir"; fi
if [ ! -d "output/kurtosis_output_dir" ] && [ "$dki" -eq 1 ]; then mkdir "output/kurtosis_output_dir"; fi
if [ ! -d "output/dsi_output_dir" ] && [ "$dsi" -eq 1 ]; then mkdir "output/dsi_output_dir"; fi

dwi=$(find . -name "*.nii.gz" -exec sh -c 'for f; do [ "$(fslval "$f" dim4)" -gt 1 ] && echo "$f"; done' sh {} +)
bval=$(find . -name "*.bval")
bvec=$(find . -name "*.bvec")

dipy_patch_rad=2
if [[ "$organism" == "mouse" ]]; then dipy_patch_rad=3; fi
if [[ "$organism" == "rat" ]]; then dipy_patch_rad=3; fi
if [[ "$organism" == "human" ]]; then dipy_patch_rad=16; fi

echo $dwi
fslmaths $dwi -Tmean -bin output/initial_dti_fullbinmask.nii.gz
dtifit --data=$dwi --mask=output/initial_dti_fullbinmask.nii.gz --bvecs=$bvec --bvals=$bval --out=output/initial_dti

t2_dat="output/initial_dti_S0.nii.gz"
brainmask="output/initial_dti_BRAINMASK.nii.gz"
headmask="output/initial_dti_HEADMASK.nii.gz"

mkdir output/temp

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Making a brain mask"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Brain mask step 1 - Bias correcting"

if [ "$organism" != "human" ]; then

    # FAST for bias correcting 
    fslmaths $t2_dat -thrP 10 output/temp/masked.nii.gz
    p95=$(fslstats output/temp/masked.nii.gz -P 95)
    fslmaths output/temp/masked.nii.gz -min ${p95//.*} ${t2_dat//.nii.gz}_clip.nii.gz
    fast -l 1 -I 20 -B -n 2 ${t2_dat//.nii.gz}_clip.nii.gz
    rm ${t2_dat//.nii.gz}_clip_seg.nii.gz ${t2_dat//.nii.gz}_clip_pve*.nii.gz ${t2_dat//.nii.gz}_clip_mixeltype.nii.gz
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Brain mask step 2 - Extracting brain"
    # Take square root and do initial thresholding at 40th percentile
    fslmaths ${t2_dat//.nii.gz}_clip_restore.nii.gz -sqrt output/temp/sqrt.nii.gz
    fslmaths output/temp/sqrt.nii.gz -thrP 40 output/temp/thresh.nii.gz
    
    # Binarize, smooth, threshold and binarize to get a mask excluding peripheral voxels which have survived thresholding, then smooth to create a cloud and multiply by original image
    fslmaths output/temp/thresh.nii.gz -bin -s 1.5 -thrP 80 -bin -s 1.5 output/temp/cloud.nii.gz
    fslmaths output/temp/cloud.nii.gz -mul ${t2_dat//.nii.gz}_clip_restore.nii.gz -thrP 20 output/temp/cloud_intermed.nii.gz
    fslmaths output/temp/cloud_intermed.nii.gz -bin -ero -ero -ero -s 0.5 -thr 0.1 -bin -dilD -dilD output/temp/cloud_mask.nii.gz
    fslmaths ${t2_dat//.nii.gz}_clip_restore.nii.gz -mas output/temp/cloud_mask.nii.gz -thrP 30 -bin -ero -ero -dilD -dilD $brainmask
    
    rm -r output/temp

else

    bet $t2_dat ${t2_dat//.nii.gz}_bet.nii.gz -m
    mv ${t2_dat//.nii.gz}_bet_mask.nii.gz $brainmask

fi

fslmaths $t2_dat -s 1 -thrP 15 -bin $headmask

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Initializing script for slice removal and MP-PCA"
/opt/conda/envs/diffusionmritoolkit/bin/python ${scripts_dir}/ODTB-1b_filt-mppca.py ${rootdir}/${d} $dwi $bval $bvec output/initial_dti $do_mppca $do_slice $dipy_patch_rad

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
        dtifit --data=$dwi --mask=$mask --bvecs=$bvec --bvals=$bval --out=output/dti_output_dir/dti
        
        #######################################
        # RUN MDT
        #######################################
        echo ""
        if [[ "$noddi" == 1 ]]; then
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Fitting NODDI to data with MDT to slice ${slice_i}"
            echo ""
            mdt-create-protocol $bvec $bval -o $protocol
            #/opt/conda/envs/diffusionmritoolkit/bin/python ${scripts_dir}/ODTB-1d_mdt.py $dwi $mask $protocol
            
            mdt-model-fit -o output/noddi_output_dir/ "NODDI" $dwi $protocol $mask
        fi
        if [[ "$dki" == 1 ]]; then
            echo ""
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Fitting DKI to data with DIPY to slice ${slice_i}"
            echo ""
            /opt/conda/envs/diffusionmritoolkit/bin/python ${scripts_dir}/ODTB-1c_dki.py output/kurtosis_output_dir/ $dwi $bval $bvec $mask
        fi
        if [[ "$dsi" == 1 ]]; then
            echo ""
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Using DSI to fit ${slice_i}"
            
            dsi_studio --action=src --source=$dwi --bval=$bval --bvec=$bvec --output=output/dsi_output_dir/dsi_out.sz --other_output=qa
            dsi_studio --action=rec --source=output/dsi_output_dir/dsi_out.sz --mask=$mask --output=output/dsi_output_dir/dsi_out.fib.gz
            dsi_studio --action=exp --source=output/dsi_output_dir/dsi_out.fib.gz  --export=qa
            /opt/conda/envs/diffusionmritoolkit/bin/python ${scripts_dir}/ODTB-1e_align.py output/dsi_output_dir/*.qa.nii.gz $mask output/dsi_output_dir/qa_align.nii.gz
            #For newer versions of dsi studio
            #dsi_studio --action=ana --source=output/dsi_output_dir/dsi_out.fib.gz  --export=qa
        fi
        echo ""
        
        if [ ! -d "output/final_output_dir" ]; then 
        mkdir "output/final_output_dir"
        cp output/dti_output_dir/dti_FA.nii.gz output/final_output_dir/FA.nii.gz
        cp output/dti_output_dir/dti_MD.nii.gz output/final_output_dir/MD.nii.gz
        cp output/dti_output_dir/dti_S0.nii.gz output/final_output_dir/S0.nii.gz
        cp output/dti_output_dir/dti_L1.nii.gz output/final_output_dir/AD.nii.gz
        fslmaths output/dti_output_dir/dti_L2.nii.gz -add output/dti_output_dir/dti_L3.nii.gz -div 2 output/final_output_dir/RD.nii.gz
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
        if [[ "$dsi" == 1 ]]; then
            cp output/dsi_output_dir/qa_align.nii.gz output/final_output_dir/QA.nii.gz
            rm -r output/dsi_output_dir
            mkdir output/dsi_output_dir
        fi
        else
        fslmaths output/final_output_dir/FA.nii.gz -add output/dti_output_dir/dti_FA.nii.gz output/final_output_dir/FA.nii.gz
        fslmaths output/final_output_dir/MD.nii.gz -add output/dti_output_dir/dti_MD.nii.gz output/final_output_dir/MD.nii.gz
        fslmaths output/final_output_dir/S0.nii.gz -add output/dti_output_dir/dti_S0.nii.gz output/final_output_dir/S0.nii.gz
        fslmaths output/final_output_dir/AD.nii.gz -add output/dti_output_dir/dti_L1.nii.gz output/final_output_dir/AD.nii.gz
        fslmaths output/dti_output_dir/dti_L2.nii.gz -add output/dti_output_dir/dti_L3.nii.gz -div 2 output/dti_output_dir/RD.nii.gz
        fslmaths output/final_output_dir/RD.nii.gz -add output/dti_output_dir/RD.nii.gz output/final_output_dir/RD.nii.gz
        rm -r output/dti_output_dir
        mkdir output/dti_output_dir
        if [[ "$noddi" == 1 ]]; then
            fslmaths output/final_output_dir/ODI.nii.gz -add output/noddi_output_dir/NODDI/ODI.nii.gz output/final_output_dir/ODI.nii.gz
            fslmaths output/final_output_dir/NDI.nii.gz -add output/noddi_output_dir/NODDI/NDI.nii.gz output/final_output_dir/NDI.nii.gz
            rm -r output/noddi_output_dir
            mkdir output/noddi_output_dir
        fi
        if [[ "$dki" == 1 ]]; then
            fslmaths output/final_output_dir/MK.nii.gz -add output/kurtosis_output_dir/MK.nii.gz output/final_output_dir/MK.nii.gz
            rm -r output/kurtosis_output_dir
            mkdir output/kurtosis_output_dir
        fi
        if [[ "$dsi" == 1 ]]; then
            fslmaths output/final_output_dir/QA.nii.gz -add output/dsi_output_dir/qa_align.nii.gz output/final_output_dir/QA.nii.gz
            rm -r output/dsi_output_dir
            mkdir output/dsi_output_dir
        fi
        fi
        
    done
    
    mv slices/filt_log.csv output/final_output_dir/filt_log.csv
    rm -r slices
    
else
    
    mask=$headmask
    protocol=protocol.prtcl
    
    #######################################
    # DTIFIT
    #######################################
    echo ""
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Fitting DTI to data with dtifit"
    echo ""
    dtifit --data=$dwi --mask=$mask --bvecs=$bvec --bvals=$bval --out=output/dti_output_dir/dti
    
    #######################################
    # RUN MDT
    #######################################
    if [[ "$noddi" == 1 ]]; then
        echo ""
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Fitting NODDI to data with MDT"
        echo ""
        mdt-create-protocol $bvec $bval -o $protocol
        mdt-model-fit -o output/noddi_output_dir/ "NODDI" $dwi $protocol $mask
    fi
    
    if [[ "$dki" == 1 ]]; then
        echo ""
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Fitting DKI to data with DIPY"
        echo ""
        /opt/conda/envs/diffusionmritoolkit/bin/python ${scripts_dir}/ODTB-1c_dki.py output/kurtosis_output_dir/ $dwi $bval $bvec $mask
        echo ""
    fi
    if [[ "$dsi" == 1 ]]; then
        echo ""
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Using DSI to fit QA"
        
        dsi_studio --action=src --source=$dwi --bval=$bval --bvec=$bvec --output=output/dsi_output_dir/${dwi//.nii.gz}.src.gz --other_output=qa
        dsi_studio --action=rec --source=output/dsi_output_dir/${dwi//.nii.gz}.src.gz --mask=$mask --output=output/dsi_output_dir/${dwi//.nii.gz}.fib.gz
        dsi_studio --action=exp --source=output/dsi_output_dir/${dwi//.nii.gz}.fib.gz  --export=qa
        /opt/conda/envs/diffusionmritoolkit/bin/python ${scripts_dir}/ODTB-1e_align.py output/dsi_output_dir/*.qa.nii.gz $mask output/dsi_output_dir/qa_align.nii.gz
        #for newer DSI studio versions
        #dsi_studio --action=ana --source=output/dsi_output_dir/${dwi//.nii.gz}.fib.gz  --export=qa
    fi
    
    mkdir "output/final_output_dir"
    cp output/dti_output_dir/dti_FA.nii.gz output/final_output_dir/FA.nii.gz
    cp output/dti_output_dir/dti_MD.nii.gz output/final_output_dir/MD.nii.gz
    cp output/dti_output_dir/dti_S0.nii.gz output/final_output_dir/S0.nii.gz
    cp output/dti_output_dir/dti_L1.nii.gz output/final_output_dir/AD.nii.gz
    fslmaths output/dti_output_dir/dti_L2.nii.gz -add output/dti_output_dir/dti_L3.nii.gz -div 2 output/final_output_dir/RD.nii.gz
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
    if [[ "$dsi" == 1 ]]; then
        cp output/dsi_output_dir/qa_align.nii.gz output/final_output_dir/QA.nii.gz
        #rm -r output/dsi_output_dir
        #mkdir output/dsi_output_dir
    fi
        
fi

mv output/final_output_dir ${d}_output
mv ${d}_output ../batch_output/${d}_output






