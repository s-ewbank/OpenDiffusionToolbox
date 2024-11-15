rootdir=${1}
vol_dir=${2}
T2_struct_vol=${3}
pre_reg_mask=${4}

cd $rootdir
cd $vol_dir

mkdir ${rootdir}/${vol_dir}/temp

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Extracting brain"
fslmaths $T2_struct_vol -sqrt ${rootdir}/${vol_dir}/temp/sqrt.nii.gz
fslmaths ${rootdir}/${vol_dir}/temp/sqrt.nii.gz -thrP 40 ${rootdir}/${vol_dir}/temp/thresh.nii.gz

com_y=$(fslstats ${rootdir}/${vol_dir}/temp/thresh.nii.gz -C | awk '{print $2}')
com_y_int=$(printf "%.0f\n" "$com_y")
fslroi ${rootdir}/${vol_dir}/temp/sqrt.nii.gz ${rootdir}/${vol_dir}/temp/sqrt_bottom.nii.gz 0 -1 0 $com_y_int 0 -1
fslroi ${rootdir}/${vol_dir}/temp/sqrt.nii.gz ${rootdir}/${vol_dir}/temp/sqrt_top.nii.gz 0 -1 $com_y_int -1 0 -1

fslmaths ${rootdir}/${vol_dir}/temp/sqrt_top.nii.gz -thrP 40 -bin -ero -ero ${rootdir}/${vol_dir}/temp/thresh_bin_top.nii.gz
fslmaths ${rootdir}/${vol_dir}/temp/sqrt_bottom.nii.gz -thrP 30 -bin -ero -ero ${rootdir}/${vol_dir}/temp/thresh_bin_bottom.nii.gz

fslmerge -y ${rootdir}/${vol_dir}/temp/thresh_bin.nii.gz ${rootdir}/${vol_dir}/temp/thresh_bin_bottom.nii.gz ${rootdir}/${vol_dir}/temp/thresh_bin_top.nii.gz

fslmaths ${rootdir}/${vol_dir}/temp/thresh_bin.nii.gz -s 0.2 -thrP 30 -bin -s 1.5 -thrP 80 -bin -s 2 ${rootdir}/${vol_dir}/temp/cloud.nii.gz

fslmaths ${rootdir}/${vol_dir}/temp/cloud.nii.gz -mul $T2_struct_vol -thrP 20 -bin -ero -ero -ero -s 0.5 -thr 0.1 -bin -dilD -dilD ${rootdir}/${vol_dir}/temp/cloud_mask.nii.gz
fslmaths $T2_struct_vol -mas ${rootdir}/${vol_dir}/temp/cloud_mask.nii.gz -thrP 30 -bin -ero -ero -dilD -dilD $pre_reg_mask
fslmaths $T2_struct_vol -mas $pre_reg_mask ${T2_struct_vol//.nii.gz}_masked.nii.gz

rm -r ${rootdir}/${vol_dir}/temp
