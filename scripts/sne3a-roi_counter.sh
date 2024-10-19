####################################################
# SET ROOTDIR - CHANGE AS NEEDED!
####################################################
rootdir="/scratch/groups/rairan/NODDI/24-08-31/batch_volumes_raw_240831/batch_output/batch_miracl_output/"
output_dir="${rootdir}/roi_analysis_output/"

sing_path=/scratch/groups/rairan/NODDI/MIRACL/miracl_latest.sif

####################################################
# NAVIGATE ROOTDIR AND SET PATHS
####################################################
cd $rootdir
subdirs=$(ls)
n_subdirs=`find . -mindepth 1 -maxdepth 1 -type d | wc -l`

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Identified $n_subdirs subdirectories in root directory for registration in $rootdir"
echo " "

mkdir $output_dir
mkdir "${output_dir}/roi_masks/"

labels_path="${rootdir}/roi_analysis_output/annotation_hemi_combined_50um_parent-level_2.nii.gz"
singularity exec $sing_path cp /code/atlases/ara/annotation/annotation_hemi_combined_50um_parent-level_2.nii.gz $labels_path

lut_path="${rootdir}/roi_analysis_output/ara_mouse_structure_graph_hemi_combined.csv"
singularity exec $sing_path cp /code/atlases/ara/ara_mouse_structure_graph_hemi_combined.csv $lut_path

####################################################
# EXECUTE ROI COUNTER
####################################################
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running script to make a bunch of little masks."
singularity exec $sing_path python /home/groups/rairan/NODDI/sne3c-roi_mask_maker.py $lut_path $labels_path $output_dir

count=1
#for d in $subdirs
#do
#    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running roi_counter script on subdir $count called $d."
#    sbatch /home/groups/rairan/NODDI/sne3b-roi_counter.sbatch $rootdir $d $lut_path $labels_path $output_dir
#    count=$(($count+1))
#done

singularity exec $sing_path python /home/groups/rairan/NODDI/sne3e-roi_out_combiner.py $output_dir
find "$output_dir" -type f ! -name '*combined*' -print0 | xargs -0 -I {} rm {}


