import os
import pandas as pd
import sys

#Input variables
output_path = sys.argv[1]

volume_types=["ODI","NDI","MD","FA"]

for vt in volume_types:
    files = [i for i in os.listdir(output_path) if (str("_"+vt+".csv") in i)]
    
    dfs=[]
    
    for f in range(len(files)):
        if f==0:
            dfs.append(pd.read_csv(output_path+"/"+files[f],index_col=0,usecols=[1,2,3,4]))
        else:
            dfs.append(pd.read_csv(output_path+"/"+files[f],index_col=0,usecols=[1,3,4]))
    final_df = pd.concat(dfs, axis=1)
    final_df.to_csv(output_path+"/combined_"+vt+".csv")