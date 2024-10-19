import nibabel as nib
import numpy as np
import os
import pandas as pd
import sys

#Input variables
lut_path = sys.argv[1]
labels_path = sys.argv[2]
output_path = sys.argv[3]

roi_lut=pd.read_csv(lut_path)
roi_img = nib.load(labels_path)

roi_data = roi_img.get_fdata()
unique_rois = np.unique(roi_data)
unique_rois = unique_rois[unique_rois != 0]

for roi in unique_rois:
    roi_mask = roi_data == roi
    try:
        roi_name = roi_lut[roi_lut.id == int(roi)].name.item()
        mask_nifti = nib.Nifti1Image(roi_mask.astype(np.int), affine=roi_img.affine)
        mask_filename = output_path+"roi_masks/"+roi_name+"_mask.nii.gz"
        nib.save(mask_nifti, mask_filename)
    except:
        continue