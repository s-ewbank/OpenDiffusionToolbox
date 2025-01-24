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
volume_types_in = sys.argv[5:]

print(volume_types_in)

print("ROI counter script")


roi_lut=pd.read_csv(lut_path)
roi_img = nib.load(labels_path)

dataset={"dir_path": dir_path}

volumes=os.listdir(dataset["dir_path"])
volume_types=[]
for vol in volume_types_in:
    if len([i for i in volumes if "_"+vol+"_reg.nii.gz" in i])>0:
        dataset[vol] = [i for i in volumes if "_"+vol+"_reg.nii.gz" in i][0]
        volume_types.append(vol)
        print(vol+" found at "+dataset[vol])
    else:
        print(vol+" not found! skipping")

# Get mask
brain_mask = nib.load(dataset["dir_path"]+"/atlas_mask.nii.gz")
brain_mask_data = brain_mask.get_fdata()

all_roi_data = roi_img.get_fdata()

if len(all_roi_data.shape)==3:
    roi_data=all_roi_data
    unique_rois = np.unique(roi_data)
    unique_rois = unique_rois[unique_rois != 0]
    
    for volume_type in volume_types:
        img = nib.load(dataset["dir_path"]+"/"+dataset[volume_type])
        img_data = img.get_fdata()
        
        roi_averages = {}
        roi_nVox = {}
        for roi in unique_rois:
            roi_mask = roi_data == roi
            masked_img = np.where(roi_mask, img_data, np.nan)
            masked_img = np.where(brain_mask_data, masked_img, np.nan)
            average_img = np.nanmean(masked_img)
            roi_nVox[roi] = np.count_nonzero(~np.isnan(masked_img))
            roi_averages[roi] = average_img
            
        data_list = []
        
        for roi, average_img in roi_averages.items():
            try:
                roi_name = roi_lut[roi_lut.id == int(roi)].name.item()
                data_list.append({"ROI_ID": int(roi), "ROI_Name": roi_name, 
                    dataset["dir_path"].split(".")[0]+"_"+volume_type: average_img})
            except:
                continue
        df = pd.DataFrame(data_list)
        df.to_csv(output_path+"/"+dataset["dir_path"].split(".")[0]+"_"+volume_type+".csv")
        
if len(all_roi_data.shape)==4:
    for volume_type in volume_types:
        img = nib.load(dataset["dir_path"]+"/"+dataset[volume_type])
        img_data = img.get_fdata()
        full_df=[]
        for atl in range(all_roi_data.shape[3]):
            roi_data=all_roi_data[:,:,:,atl]
            unique_rois = np.unique(roi_data)
            unique_rois = unique_rois[unique_rois != 0]
            
            roi_averages = {}
            roi_nVox = {}
            for roi in unique_rois:
                roi_mask = roi_data == roi
                masked_img = np.where(roi_mask, img_data, np.nan)
                masked_img = np.where(brain_mask_data, masked_img, np.nan)
                average_img = np.nanmean(masked_img)
                roi_nVox[roi] = np.count_nonzero(~np.isnan(masked_img))
                roi_averages[roi] = average_img
                
            data_list = []
            
            for roi, average_img in roi_averages.items():
                try:
                    roi_name = roi_lut[roi_lut.id == int(roi)].name.item()
                    data_list.append({"ROI_ID": int(roi), "ROI_Name": roi_name, 
                        dataset["dir_path"].split(".")[0]+"_"+volume_type: average_img})
                except:
                    continue
            full_df.append(pd.DataFrame(data_list))
        df=pd.concat(full_df,axis=0)
        df.to_csv(output_path+"/"+dataset["dir_path"].split(".")[0]+"_"+volume_type+".csv")