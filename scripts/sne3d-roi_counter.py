import nibabel as nib
import numpy as np
import os
import pandas as pd
import sys

#Input variables
dir_path = sys.argv[1]
lut_path = sys.argv[2]
labels_path = sys.argv[3]
output_path = sys.argv[4]

roi_lut=pd.read_csv(lut_path)
dataset={"dir_path": dir_path}
roi_img = nib.load(labels_path)

volumes=os.listdir(dataset["dir_path"])
dataset["ODI"] = [i for i in volumes if "_ODI_reg.nii.gz" in i][0]
dataset["NDI"] = [i for i in volumes if "_NDI_reg.nii.gz" in i][0]
dataset["MD"] = [i for i in volumes if "_MD_reg.nii.gz" in i][0]
dataset["FA"] = [i for i in volumes if "_FA_reg.nii.gz" in i][0]
volume_types=["ODI","NDI","MD","FA"]

# Get mask
brain_mask = nib.load(dataset["dir_path"]+"/miracl_brain_mask.nii.gz")
brain_mask_data = brain_mask.get_fdata()

roi_data = roi_img.get_fdata()
unique_rois = np.unique(roi_data)
unique_rois = unique_rois[unique_rois != 0]

count=0
for volume_type in volume_types:
    img = nib.load(dataset["dir_path"]+"/"+dataset[volume_type])
    img_data = img.get_fdata()
    
    roi_averages = {}
    roi_nVox = {}
    for roi in unique_rois:
        roi_mask = roi_data == roi
        masked_img = np.where(roi_mask, img_data, np.nan)
        masked_img = np.where(brain_mask_data, masked_img, np.nan)
        average_img = np.nansum(masked_img)
        roi_nVox[roi] = np.count_nonzero(~np.isnan(masked_img))
        roi_averages[roi] = average_img
        
    count=0
    data_list = []
    
    for roi, average_img in roi_averages.items():
        try:
            roi_name = roi_lut[roi_lut.id == int(roi)].name.item()
            data_list.append({"ROI_ID": int(roi), "ROI_Name": roi_name, 
                dataset["dir_path"].split(".")[0]+"_nVox": roi_nVox[roi], 
                dataset["dir_path"].split(".")[0]+"_"+volume_type: average_img})
        except:
            continue
    df = pd.DataFrame(data_list)
    df.to_csv(output_path+"/"+dataset["dir_path"].split(".")[0]+"_"+volume_type+".csv")
