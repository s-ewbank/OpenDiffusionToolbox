import os
import sys
import nibabel as nib
import nilearn
from nilearn import plotting, decoding, image
import numpy as np

path = sys.argv[1]
volume_type = sys.argv[2]
mask = sys.argv[3]
organism = sys.argv[4]

if organism=="mouse":
    cut_coords=np.arange(-13,-1,1.5)
elif organism=="rat":
    cut_coords=np.arange(-13,6,2.5)
elif organism=="human":
    cut_coords=np.arange(-60,60,20)


volumes = [path+"/"+volume_type+"/"+i for i in os.listdir(path+"/"+volume_type)]

### New approach
print("Making average NIFTI volume")
vol0 = nib.load(volumes[0])
vol0_data = vol0.get_fdata()
n_slices = vol0_data.shape[2]
sum_data = np.zeros_like(vol0_data)

for volume in volumes:
    img = nib.load(volume)
    data = img.get_fdata()
    sum_data += data
    
avg_data = sum_data / len(volumes)
avg_img = nib.Nifti1Image(avg_data, affine=vol0.affine)
output_name=str(path+"avg_"+volume_type+".nii.gz")
nib.save(avg_img, output_name)

# OLD APPRAOCH
# print("Making average NIFTI volume")
# avg_img=nilearn.image.mean_img(volumes)
# output_name=str(path+"avg_"+volume_type+".nii.gz")
# nib.save(avg_img, output_name)

print("Masking average NIFTI")
mask=nilearn.image.load_img(mask)
avg_img_masked=nilearn.image.math_img("mask*avg",avg=avg_img,mask=mask)
masked_output_name=str(path+"avg_"+volume_type+"_masked.nii.gz")
nib.save(avg_img_masked, masked_output_name)


nilearn.plotting.plot_img(avg_img,title="Average "+volume_type,black_bg=True,cmap="gray",
    draw_cross=False,annotate=False,output_file=path+"avg_"+volume_type+".png",
    display_mode="y",cut_coords=cut_coords)
    

nilearn.plotting.plot_img(avg_img_masked,title="Average "+volume_type,black_bg=True,cmap="gray",
    draw_cross=False,annotate=False,output_file=path+"avg_"+volume_type+"_masked.png",
    display_mode="y",cut_coords=cut_coords)