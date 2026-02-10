import pandas as pd
import numpy as np
import sys

groups_df_path = sys.argv[1]
save_path = sys.argv[2]

df = pd.read_csv(groups_df_path,header=0)
n_subjs = df.shape[0]
groups = sorted(list(df.group.unique()))
corr_data_name = [col for col in df.columns if col.startswith("data_")][0]

save_stem = str("GROUP-" + groups[0] + "_CORR-" + corr_data_name.split("data_")[1])

# Making design mat file
design_mat = df[corr_data_name].values.reshape(-1, 1)
design_mat = (design_mat - np.mean(design_mat)) / np.std(design_mat)

np.savetxt(save_path+save_stem+"_design.mat", design_mat, fmt='%0.6f')

#Making design con file
design_con = np.array([[1], [-1]])
design_con = design_con.astype(int)
np.savetxt(save_path+save_stem+"_design.con", design_con, fmt='%d')

print(save_stem)