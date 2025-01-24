import os
import pandas as pd
import sys

data_dir = sys.argv[1]

stat_files = [i for i in os.listdir(data_dir) if ".stat." in i]

all_data=[]

for f in stat_files:
    df=pd.read_csv(data_dir+"/"+f,sep='\t',header=None,index_col=0)
    df.columns=[f]
    all_data.append(df)

all_data_df = pd.concat(all_data,axis=1)

all_data_df.to_csv(data_dir+"/tract-stats.csv")