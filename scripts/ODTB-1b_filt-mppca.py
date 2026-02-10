import numpy as np
import nibabel as nib
from dipy.denoise.localpca import mppca
from datetime import datetime
import sys
import pandas as pd
import os

path = sys.argv[1]
raw_dwi_path = sys.argv[2]
bval_path = sys.argv[3]
bvec_path = sys.argv[4]
dtifit_prefix = sys.argv[5]
do_mppca = int(sys.argv[6])
do_slices = int(sys.argv[7])
patch_rad = int(sys.argv[8])

S0=nib.load(dtifit_prefix+'_S0.nii.gz')
affine=S0.affine
S0=S0.get_fdata()

L1=nib.load(dtifit_prefix+'_L1.nii.gz').get_fdata()
L2=nib.load(dtifit_prefix+'_L2.nii.gz').get_fdata()
L3=nib.load(dtifit_prefix+'_L3.nii.gz').get_fdata()

V1=nib.load(dtifit_prefix+'_V1.nii.gz').get_fdata()
V2=nib.load(dtifit_prefix+'_V2.nii.gz').get_fdata()
V3=nib.load(dtifit_prefix+'_V3.nii.gz').get_fdata()

brainmask=nib.load(dtifit_prefix+'_BRAINMASK.nii.gz').get_fdata()
headmask=nib.load(dtifit_prefix+'_HEADMASK.nii.gz').get_fdata()

bvec_orig=np.loadtxt(bvec_path)
bval=np.loadtxt(bval_path)
b0=np.sum(bval<100)

try:
    bvec=bvec_orig[bval>100]
except IndexError:
    bvec_orig=bvec_orig.T
    bvec=bvec_orig[bval>100]

bval=bval[bval>100]

def calculate_signal(bvec_in, bval_in, D, S0_in):
    b = np.array(bvec_in)
    b = b / np.linalg.norm(b) 
    b_squared = np.dot(b.T, D @ b)
    return S0_in * np.exp(-bval_in * b_squared)

def mean_squared_error(arr1,arr2):
    return np.mean(np.square(arr1-arr2))

if do_slices==1:
    os.makedirs(path+"/slices/",exist_ok=True)
    tensorperfect_dwi=np.zeros([S0.shape[0],S0.shape[1],S0.shape[2],len(bval)])
    
    n_voxels=np.sum(headmask)
    voxels_completed=0
    
    print(datetime.now().strftime("[%Y-%m-%d %H:%M] ") + f'Starting reconstruction of tensor-perfect raw diffusion image - 0% of {n_voxels} voxels completed')
    for x in range(S0.shape[0]):
        for y in range(S0.shape[1]):
            for z in range(S0.shape[2]):
                if headmask[x,y,z]!=0:
                    S0_vx = S0[x,y,z]
                    eigvec_vx = np.array([V1[x,y,z,:],V2[x,y,z,:],V3[x,y,z,:]]).T
                    eigval_vx = np.array([L1[x,y,z],L2[x,y,z],L3[x,y,z]])
                    D = eigvec_vx @ np.diag(eigval_vx) @ eigvec_vx.T
                    for dxn in range(len(bval)):
                        bvec_dxn = bvec[dxn,:]
                        bval_dxn = bval[dxn]
                        signal = calculate_signal(bvec_dxn, bval_dxn, D, S0_vx)
                        tensorperfect_dwi[x,y,z,dxn] = signal
                    voxels_completed+=1
                    if (voxels_completed/n_voxels)%0.25==0:
                        print(datetime.now().strftime("[%Y-%m-%d %H:%M] ") + f'{100*voxels_completed/n_voxels}% of voxels completed')
    
    for b in range(b0):
        tensorperfect_dwi=np.concatenate((S0[..., np.newaxis], tensorperfect_dwi), axis=3)
    tensorperfect_dwi_img=nib.Nifti1Image(tensorperfect_dwi, affine)
    nib.save(tensorperfect_dwi_img,path+'/tensorperfect_dwi.nii.gz')
    
    raw_dwi_img = nib.load(raw_dwi_path)
    raw_dwi = raw_dwi_img.get_fdata()
    
    dwi_censored=np.zeros_like(raw_dwi)
    censor_mask=np.zeros_like(raw_dwi)
    bval_orig=np.loadtxt(bval_path)
    bval=(bval_orig/500).round()*500
    bval_un=sorted(np.unique(bval))
    bval_masks={}
    for bval_i in bval_un:
        bval_masks[bval_i]=np.argwhere(bval==bval_i).T[0]
    
    slices_dropped=0
    slice_keep_log={}

    bbox_x=np.sum(brainmask,axis=(1,2))>0
    bbox_x=np.where(bbox_x)[0][[0, -1]]
    bbox_y=np.sum(brainmask,axis=(0,2))>0
    bbox_y=np.where(bbox_y)[0][[0, -1]]

    raw_dwi_bounded=raw_dwi[bbox_x[0]:bbox_x[1],bbox_y[0]:bbox_y[1],:,:]

    for slice_i in range(tensorperfect_dwi.shape[2]):
        keep_slice=[]
        for b in bval_un:
            bval_m=bval_masks[b]
            if b==0:
                for volume in bval_m:
                    dwi_censored[:,:,slice_i,volume]=raw_dwi[:,:,slice_i,volume]
                    keep_slice.append(True)
            else:
                raw_dwi_b=raw_dwi_bounded[:,:,:,bval_m]
                raw_dwi_b_mean=np.median(raw_dwi_b,axis=3)
                for volume in bval_m:
                    bads=[]
                    #mse = mean_squared_error(raw_dwi_bounded[:,:,slice_i,volume],raw_dwi_b_mean[:,:,slice_i])
                    #if np.sqrt(mse)>0.25*np.std(raw_dwi_bounded[:,:,slice_i,:]):
                    err=np.sqrt((raw_dwi_bounded[:,:,slice_i,volume]-raw_dwi_b_mean[:,:,slice_i])**2)
                    fr_bad = np.mean(err>(0.5*np.std(raw_dwi_bounded[:,:,slice_i,:],axis=2)))
                    if fr_bad>0.5:
                        keep_slice.append(False)
                        dwi_censored[:,:,slice_i,volume]=tensorperfect_dwi[:,:,slice_i,volume]
                    else:
                        keep_slice.append(True)
                        dwi_censored[:,:,slice_i,volume]=raw_dwi[:,:,slice_i,volume]

        print(f"Slice {slice_i} - {np.sum(np.array(keep_slice)==False)} bads")
        slice_keep_log[slice_i]=keep_slice
        np.savetxt(path+f'/slices/bval_{slice_i}.bval', bval_orig[keep_slice], fmt='%.2f', newline=' ')
        np.savetxt(path+f'/slices/bvec_{slice_i}.bvec', bvec_orig[keep_slice,:], fmt='%.6f', delimiter=' ', newline='\n')

    log_df = pd.DataFrame(slice_keep_log)
    log_df = log_df.astype(int)
    log_df.to_csv(path+'/slices/filt_log.csv')
    
    dwi_censored_img=nib.Nifti1Image(dwi_censored, affine)
    nib.save(dwi_censored_img,path+'/raw_dwi_censored.nii.gz')
    
    if do_mppca==1:
        print(datetime.now().strftime("[%Y-%m-%d %H:%M] ") + f'Doing MP-PCA')
        dwi_censored_mppca = mppca(dwi_censored, patch_radius=[patch_rad, patch_rad, patch_rad],mask=headmask)
        dwi_censored_mppca_img=nib.Nifti1Image(dwi_censored_mppca, affine)
        nib.save(dwi_censored_mppca_img,path+'/raw_dwi_censored_mppca.nii.gz')
        for slice_i in range(tensorperfect_dwi.shape[2]):
            keep_slice=slice_keep_log[slice_i]
            slice_dat = dwi_censored_mppca[:,:,:,keep_slice]
            slice_img = nib.Nifti1Image(slice_dat, affine)
            nib.save(slice_img, path+f'/slices/dwi_censored_mppca_{slice_i}.nii.gz')
            mask_i = np.zeros_like(headmask)
            mask_i[:,:,slice_i] = headmask[:,:,slice_i]
            mask_img = nib.Nifti1Image(mask_i, affine)
            nib.save(mask_img, path+f'/slices/mask_{slice_i}.nii.gz')
    
    else:
        for slice_i in range(tensorperfect_dwi.shape[2]):
            keep_slice=slice_keep_log[slice_i]
            slice_dat = dwi_censored[:,:,:,keep_slice]
            slice_img = nib.Nifti1Image(slice_dat, affine)
            nib.save(slice_img, path+f'/slices/dwi_censored_{slice_i}.nii.gz')

elif do_mppca==1:
    raw_dwi_img = nib.load(raw_dwi_path)
    raw_dwi = raw_dwi_img.get_fdata()
    print(datetime.now().strftime("[%Y-%m-%d %H:%M] ") + f'Doing MP-PCA')
    raw_dwi_mppca = mppca(raw_dwi, patch_radius=[patch_rad, patch_rad, patch_rad], mask=headmask)
    raw_dwi_mppca_img=nib.Nifti1Image(raw_dwi_mppca, affine)
    nib.save(raw_dwi_mppca_img,path+'/raw_dwi_mppca.nii.gz')
    print(datetime.now().strftime("[%Y-%m-%d %H:%M] ") + f'Done with MP-PCA')



