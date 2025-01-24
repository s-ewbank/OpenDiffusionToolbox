import os
import pandas as pd
import sys

#Input variables
output_path = sys.argv[1]
volume_types = sys.argv[2:]
print(volume_types)

for vt in volume_types:
    files = [i for i in os.listdir(output_path) if (str("_"+vt+".csv") in i)]
    
    if len(files)>0:
        dfs=[]
        
        for f in range(len(files)):
            if f==0:
                dfs.append(pd.read_csv(output_path+"/"+files[f],index_col=0,usecols=[1,2,3]))
            else:
                dfs.append(pd.read_csv(output_path+"/"+files[f],index_col=0,usecols=[1,3]))
        final_df = pd.concat(dfs, axis=1)
        final_df.to_csv(output_path+"/combined_"+vt+".csv")
    else:
        print("No files found for volume type "+vt+", continuing.")