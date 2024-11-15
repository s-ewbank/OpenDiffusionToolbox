cd /scratch/users/snewbank/MKET2/data/batch_DSI_output
dsi_sing_path=/home/groups/rairan/NODDI/DSI/dsistudio_unlhcc.sif
subdirs=$(ls)
for d in $subdirs; do
    cd $d
    bs=$(basename $d)
    tt=$(find . -type f -name "*_DSI.tt.gz")
    fib=$(find . -type f -name "*_DSI.fib.gz")
    atlas=$(find . -type f -name "atlas_warp_abbrev.nii.gz")
    ODI=/scratch/users/snewbank/MKET2/data/batch_output/${bs//_DSI}_output/ODI.nii.gz
    NDI=/scratch/users/snewbank/MKET2/data/batch_output/${bs//_DSI}_output/NDI.nii.gz

    singularity exec $dsi_sing_path dsi_studio --action=ana \
        --source=$fib \
        --tract=$tt \
        --connectivity=atlas_warp_abbrev.nii.gz \
        --other_slices=${ODI},${NDI} \
        --connectivity_value=ODI,NDI \
        --connectivity_output=connectogram
    cd ..
done