import nibabel as nib
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from time import time

from dipy.denoise.localpca import mppca
from dipy.io.image import load_nifti
from dipy.core.gradients import gradient_table
from dipy.io.gradients import read_bvals_bvecs
from dipy.segment.mask import median_otsu
import dipy.reconst.dki as dki
import sys
from scipy.ndimage import gaussian_filter

path = sys.argv[1]

def bad_slices_out(orig_image_path,
                   proc_image_path,
                   bvec_path,
                   bval_path, 
                   type="save_chopped_img",
                   lo_q=0.1,
                   hi_q=0.9):
    
    bval=np.loadtxt(bval_path)
    bval_clip=np.round(bval/500)*500
    
    b0=np.argwhere(bval_clip==0).T[0]
    b1000=np.argwhere(bval_clip==1000).T[0]
    b2500=np.argwhere(bval_clip==2500).T[0]
    
    bvec=np.loadtxt(bvec_path)
    
    img = nib.load(orig_image_path)
    proc_img = nib.load(proc_image_path)
    data = img.get_fdata()
    proc_data = proc_img.get_fdata()
    
    x_dim, y_dim, z_dim, n_vol = data.shape

    chopped_img_data=proc_data.copy()
    
    os.makedirs(path+"/slices/",exist_ok=True)
    
    slices_dropped=0
    slice_keep_log={}
    for slice_i in range(z_dim):
        new_img_data=np.zeros_like(data)
        keep_slice=[]
        for bval_i in (b0,b1000,b2500):
            hi_q_img=np.quantile(data[:,:,slice_i,bval_i],hi_q,axis=2)
            lo_q_img=np.quantile(data[:,:,slice_i,bval_i],lo_q,axis=2)
            for i in bval_i:
                slice=data[:,:,slice_i,i]
                slice_hi=slice>hi_q_img
                slice_lo=slice<lo_q_img
                slice_oddvox=slice_hi+slice_lo
                if np.mean(slice_oddvox)<0.25:
                    new_img_data[:,:,slice_i,i]=proc_data[:,:,slice_i,i]
                    keep_slice.append(True)
                else:
                    slices_dropped+=1
                    keep_slice.append(False)
                    chopped_img_data[:,:,slice_i,i]=np.zeros_like(chopped_img_data[:,:,slice_i,i])
                    #print(f"Volume {i}, slice {slice_i} - {np.mean(slice_oddvox)}")

        if type=="save_slices":
            new_img_data = new_img_data[:,:,:,keep_slice]
            np.savetxt(path+f'/slices/bval_{slice_i}.bval', bval[keep_slice], fmt='%.2f', newline=' ')
            np.savetxt(path+f'/slices/bvec_{slice_i}.bvec', bvec[keep_slice,:], fmt='%.6f', delimiter=' ', newline='\n')
            new_img = nib.Nifti1Image(new_img_data, img.affine)
            nib.save(new_img, path+f'/slices/filt_{slice_i}.nii.gz')
            slice_keep_log[slice_i]=keep_slice
        
    prop_slices_dropped=slices_dropped/(n_vol*z_dim)
    print(f"Dropped {prop_slices_dropped} of slices")

    if type=="save_chopped_img":
        chopped_new_img = nib.Nifti1Image(chopped_img_data, img.affine)
        nib.save(chopped_new_img, path+'/filt.nii.gz')
    elif type=="save_slices":
        log_df = pd.DataFrame(slice_keep_log)
        log_df = log_df.astype(int)
        log_df.to_csv(path+'/slices/filt_log.csv')

orig_image_path=[path+"/"+i for i in os.listdir(path) if ((i.endswith(".nii.gz")) and ("1_DTI" in i))][0]
bvec_path=[path+"/"+i for i in os.listdir(path) if (i.endswith(".bvec"))][0]
bval_path=[path+"/"+i for i in os.listdir(path) if (i.endswith(".bval"))][0]

print(f"Bval path is {bval_path}")
print(f"Bvec path is {bvec_path}")

bad_slices_out(orig_image_path,
               orig_image_path, #proc_image
               bvec_path,
               bval_path,
               type="save_chopped_img")

print("MP-PCA - Loading nifti")
data, affine, img = load_nifti(path+"/filt.nii.gz", return_img=True)

print("MP-PCA - Doing mppca")
t = time()
denoised_arr = mppca(data, patch_radius=2)
print("Time taken for local MP-PCA ", -t + time())
nib.save(nib.Nifti1Image(denoised_arr, affine), path+"/filt_mppca.nii.gz")

bad_slices_out(orig_image_path,
               path+"/filt_mppca.nii.gz",
               bvec_path,
               bval_path,
               type="save_slices")
               
               

mask_img = nib.load(orig_image_path)
mask_data = mask_img.get_fdata()
x_dim, y_dim, z_dim, n_vol = mask_data.shape
mask_data_sum=np.sum(mask_data[:,:,:,0:10],axis=3)
blurred_mask_data = gaussian_filter(mask_data_sum, sigma=4)
thr=np.quantile(blurred_mask_data,0.6)
blurred_mask_data=np.array(blurred_mask_data>thr,dtype=np.int16)

for i in range(z_dim):
    mask_i=np.zeros_like(blurred_mask_data)
    mask_i[:,:,i]=blurred_mask_data[:,:,i]
    nib.save(nib.Nifti1Image(mask_i, mask_img.affine,dtype=np.int16), path+f"/slices/mask_{i}.nii.gz")
