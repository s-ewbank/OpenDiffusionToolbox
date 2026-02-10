import os
import pandas as pd
import numpy as np
import sys

#Input variables
root_dir = sys.argv[1]

tract_path = root_dir + "/batch_prob_tract_output/"
files = [i for i in os.listdir(tract_path) if os.path.isdir(tract_path+i)]
files = sorted([file_i for file_i in files if "nodif_brain_mask.nii.gz" in os.listdir(tract_path+file_i)])

for f, file_i in enumerate(files):
    if f==0:
        roi_list = [i.split("probtrackx_")[1] for i in os.listdir(tract_path+file_i) if ((os.path.isdir(tract_path+file_i+"/"+i)) and ("probtrackx_" in i))]
        dfs = {roi: {} for roi in roi_list}
    for roi in roi_list:
        roi_path = tract_path + file_i + "/probtrackx_" + roi + "/waytotal"
        with open(roi_path, 'r') as file:
            waytotal_num = file.read().strip()
            dfs[roi][file_i] = int(waytotal_num)
    
final_df = pd.DataFrame(dfs)
final_df.to_csv(tract_path + "special_roi_waytotals.csv")