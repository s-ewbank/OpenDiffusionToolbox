tbss_dir=${1}

template=/Users/snewbank/Documents/PhD/Data/NODDI/TEMPLATES_100um/avg_FA_masked.nii.gz

cd $tbss_dir

mkdir stats
mkdir FA

mv *.nii.gz FA
cp $template stats/mean_FA.nii.gz
fslmaths stats/mean_FA.nii.gz -bin stats/mean_FA_mask.nii.gz
fslmerge -t stats/all_FA `imglob FA/*.nii.gz`

tbss_skeleton -i stats/mean_FA.nii.gz -o stats/mean_FA_skeleton.nii.gz
tbss_4_prestats 0.22