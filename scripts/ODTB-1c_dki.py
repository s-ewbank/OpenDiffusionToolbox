import numpy as np
import nibabel as nib

from dipy.core.gradients import gradient_table
from dipy.io.gradients import read_bvals_bvecs
from dipy.io.image import load_nifti
import dipy.reconst.msdki as msdki
from dipy.core.gradients import gradient_table
import sys

path=sys.argv[1]
raw_dwi_path=sys.argv[2]
bval_path=sys.argv[3]
bvec_path=sys.argv[4]
mask_path=sys.argv[5]

bvec_arr=np.loadtxt(bvec_path)
bval_arr=np.loadtxt(bval_path)
bval_arr[bval_arr<100]=0
gtab=gradient_table(bval_arr,bvec_arr)

dwi, affine = load_nifti(raw_dwi_path)
mask, maffine = load_nifti(mask_path)

try:
    msdki_model = msdki.MeanDiffusionKurtosisModel(gtab)
    msdki_fit = msdki_model.fit(dwi, mask=mask)
    msk = msdki_fit.msk
    
except:
    msk = np.zeros_like(dwi[:,:,:,0])

mk_img=nib.Nifti1Image(msk, affine)

nib.save(mk_img, path+'/MK.nii.gz')
