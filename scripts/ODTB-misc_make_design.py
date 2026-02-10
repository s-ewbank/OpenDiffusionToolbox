import pandas as pd
import numpy as np
import sys

groups_df_path = sys.argv[1]
save_path = sys.argv[2]

df = pd.read_csv(groups_df_path,header=0)
groups = sorted(list(df.group.unique()))
n_groups = len(groups)
n_subjs = df.shape[0]

save_stem = str("GROUP1-" + groups[0] + "_VS-GROUP2-" + groups[1])

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

print(save_stem)
