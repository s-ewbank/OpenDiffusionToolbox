import pandas as pd
import numpy as np
import sys
import nibabel as nib

groups_df_path = sys.argv[1]
merged_save_path = sys.argv[2]
data_path = sys.argv[3]
vt = sys.argv[4]

save_path = merged_save_path+vt+"_"

df = pd.read_csv(groups_df_path,header=0)
df = df.dropna(subset=['subject', 'filename', 'group', 'timepoint'])
df['filename'] = df['filename'].str.strip()
df['subject'] = df['subject'].astype(str).str.strip()
df['group'] = df['group'].astype(str).str.strip()

groups = sorted(list(df.group.unique()))
n_groups = len(groups)
subjs = df.subject.unique()
df_new = []

for subj in subjs:
    try:
        group = df.loc[(df["timepoint"]==0) & (df["subject"]==subj)].group.values[0]
        fname_0 = df.loc[(df["timepoint"]==0) & (df["subject"]==subj)].filename.values[0]
        fname_1 = df.loc[(df["timepoint"]==1) & (df["subject"]==subj)].filename.values[0]
        df_new.append({"subject": subj, "group": group, "filename_0":fname_0, "filename_1": fname_1})
    except:
        continue

df = pd.DataFrame(df_new)
n_subjs = df.shape[0]

save_stem = str("DELTA_GROUP1-" + groups[0] + "_VS-GROUP2-" + groups[1])

# Making design mat file
design_mat = np.zeros([n_subjs,n_groups])

for g, group in enumerate(groups):
    design_mat[:,g]=df.group==group

design_mat = design_mat.astype(int)
np.savetxt(save_path+save_stem+"_design.mat", design_mat, fmt='%d')

#Making design con file
design_con = np.ones(n_groups)
design_con = design_con * -1
design_con_2 = np.eye(n_groups)
design_con_2 = design_con_2 * 2
design_con = design_con + design_con_2
design_con = design_con.astype(int)
np.savetxt(save_path+save_stem+"_design.con", design_con, fmt='%d')

merged_data = []


subjs = df.subject.unique()
got_header = False

for subj in subjs:
    f0 = df.loc[(df["subject"]==subj)].filename_0.values[0]
    f1 = df.loc[(df["subject"]==subj)].filename_1.values[0]
    if not got_header:
        vol0 = nib.load(data_path+"/"+f0+"_output_reg_output/"+f0+"_output_"+vt+"_reg.nii.gz")
        affine = vol0.affine
        header = vol0.header
    vol0 = nib.load(data_path+"/"+f0+"_output_reg_output/"+f0+"_output_"+vt+"_reg.nii.gz").get_fdata()
    vol1 = nib.load(data_path+"/"+f1+"_output_reg_output/"+f1+"_output_"+vt+"_reg.nii.gz").get_fdata()
    delta = vol1-vol0
    merged_data.append(delta)

    
merged_data = np.stack(merged_data, axis=-1)
merged_data = nib.Nifti1Image(merged_data, affine, header=header)
nib.save(merged_data, merged_save_path+"merged_"+vt+"_"+save_stem+".nii.gz")

print(save_stem)