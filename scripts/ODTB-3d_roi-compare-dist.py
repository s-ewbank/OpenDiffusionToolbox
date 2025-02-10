√###############################################
# IMPORTS
###############################################
import nilearn
import sys
from nilearn import plotting, decoding, image
from nilearn import regions
import nibabel as nib
import pandas as pd
import numpy as np
import os
from scipy.stats import ttest_ind, kstest
from scipy.ndimage import label, generate_binary_structure
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import ListedColormap
rcParams['font.family'] = 'sans-serif'
plt.rc('text', usetex=False)
rcParams['font.sans-serif'] = ['Arial']
import numpy as np
from statsmodels.stats.multitest import multipletests
import argparse


###############################################
# INITIALIZE FILE, DIRECTORY NAMES
###############################################
parser = argparse.ArgumentParser(description="Script for generating stat maps.")

parser.add_argument("--rootdir", required=True, help="root directory path")
parser.add_argument("--mask", required=True, help="path to mask")
parser.add_argument("--outdir", required=True, help="output directory path")
parser.add_argument("--groups", required=True, help="path to groups file")
parser.add_argument("--organism", required=True, help="organism")
parser.add_argument("--atlas_path", required=True, help="path to atlas directory")
parser.add_argument("--vt", required=True, help="name of volume type (e.g. FA, MD)")
parser.add_argument("--design", required=True, help="statistical design; baseline, delta, or longitudinal")

args = parser.parse_args()
rootdir = args.rootdir
mask_path = args.mask
output_dir = args.outdir
groups_file = args.groups
organism = args.organism
atlas_path = args.atlas_path
vt = args.vt
design = args.design

smooth=True

if organism=="mouse":
    cut_coords=np.arange(-5,5,1.5)
    fwhm=0.1
    atlas_path=atlas_path+"/mouse/atlas_levels/ATLAS_LVL7_100um.nii.gz"
elif organism=="rat":
    cut_coords=np.arange(-13,6,2.5)
    fwhm=0.5
    atlas_path=atlas_path+"/rat/WHS_SD_rat_BRAIN_ATLAS.nii.gz"
elif organism=="human":
    cut_coords=np.arange(-60,60,20)
    fwhm=5
    atlas_path=atlas_path+"/human/HarvardOxford-sub-maxprob-thr50-1mm.nii.gz"

groups_df=pd.read_csv(groups_file)
groups=sorted(groups_df["group"].unique())

timepoints=sorted(groups_df["timepoint"].unique())

no_timepts_dict={}
for group in groups:
    group_files=groups_df[groups_df["group"]==group]["filename"].to_list()
    group_files=[i+"_output_reg_output" for i in group_files]
    no_timepts_dict[group]=group_files

timepts_dict={}
for group in groups:
    sub_group_df=groups_df[groups_df["group"]==group]
    subjs=sub_group_df["subject"].unique()
    group_timepts=[]
    for subj in subjs:
        subj_timepts=[]
        subj_df=sub_group_df[sub_group_df["subject"]==subj]
        for timept in timepoints:
            file=subj_df[subj_df["timepoint"]==timept]["filename"].to_list()
            if len(file)==0:
                file=None
                print("Warning - timepoint "+str(timept)+" not fount for subj "+str(subj))
            else:
                file=file[0]+"_output_reg_output"
            subj_timepts.append(file)
        group_timepts.append(subj_timepts)
    timepts_dict[group]=group_timepts


def compare_roi_dist(roi,groups,vt,atlas_path):
    atlas_im = nib.load(atlas_path)
    atlas_data = atlas_im.get_fdata()
    groups_data=[]
    groups_data_point=[]
    groups=sorted(list(timepts_dict.keys()))

    for group in groups:
        print("Group "+str(group)+" contains:")
        [print(i) for i in timepts_dict[group]]
        all_data_dist=[]
        all_data_point=[]
        image_pairs=timepts_dict[group]
        for subj in range(len(image_pairs)):
            subj_id=image_pairs[subj][0].split("_")[0]
            try:
                subj_pre_dir=image_pairs[subj][0]
                subj_pre_vol=[rootdir+"/"+subj_pre_dir+"/"+i for i in os.listdir(rootdir+"/"+subj_pre_dir) if "_"+vt+"_reg.nii.gz" in i][0]
            except:
                print("Subj " + subj_id + " had some issue, so skipping.")
                continue
            if smooth==True:
                subj_pre_vol=nib.load(subj_pre_vol)
                subj_pre_vol=nilearn.image.smooth_img(subj_pre_vol, fwhm)
            else:
                subj_pre_vol=nib.load(subj_pre_vol)
            
            im_data=subj_pre_vol.get_fdata()
            all_data_dist.append(im_data[atlas_data==roi])
            all_data_point.append(np.mean(im_data[atlas_data==roi]))
        all_data_dist=np.array(all_data_dist).flatten()
        groups_data.append(all_data_dist)
        groups_data_point.append(all_data_point)
    grp1_greater=np.mean(groups_data_point[0])<np.mean(groups_data_point[1])
    if grp1_greater:
        dxn=1
    else:
        dxn=-1
    result = kstest(groups_data[0],groups_data[1])
    return result, dxn

###############################################
# GETTING SUBJ IDS AND PLOTTING DIFFS AND T SCORES
###############################################
if design=="baseline":

    cmap="PiYG_r"
    #cmap="BrBG_r"
    
    atlas_im = nib.load(atlas_path)
    atlas_data = atlas_im.get_fdata()
    result_atlas=np.zeros_like(atlas_data)
    pval_atlas=np.zeros_like(atlas_data)
    for roi in np.unique(atlas_data):
        if roi!=0:
            result, dxn = compare_roi_dist(roi,groups,vt,atlas_path)
            result_atlas[atlas_data==roi]=result.statistic*dxn
            pval_atlas[atlas_data==roi]=result.pvalue
    
    result_atlas_im = nib.Nifti1Image(result_atlas, atlas_im.affine, atlas_im.header)
    pval_atlas_im = nib.Nifti1Image(pval_atlas, atlas_im.affine, atlas_im.header)
    
    outline_im = make_outline_im(atlas_path)
    
    fig, ax = plt.subplots(1,1,figsize=(8,2),dpi=500)
    
    clear_black_cmap = ListedColormap([(0, 0, 0, 0),(0, 0, 0, 0.5)], name="clear_black")

    plotting.plot_img(result_atlas_im,black_bg=False,cmap=cmap,title=groups[0]+" vs "+groups[1] + " " + vt,
                      draw_cross=False,annotate=False,display_mode="y",cut_coords=cut_coords,vmin=-1,vmax=1,axes=ax)
    
    plotting.plot_img(outline_im,black_bg=False,cmap=clear_black_cmap,title=groups[0]+" vs "+groups[1] + " " + vt,
                      draw_cross=False,annotate=False,display_mode="y",cut_coords=cut_coords,vmin=0,vmax=1,axes=ax)
                      
    plt.savefig(output_dir+"/roi-ks-result-"+groups[0]+" vs "+groups[1] + " " + vt + ".png", dpi=500)
    
    result_atlas_im.to_filename(output_dir + "/roi-ks-STAT_"+groups[0]+"_vs_"+groups[1] + "-" + vt + ".nii.gz")
    result_atlas_im.to_filename(output_dir + "/roi-ks-PVAL_"+groups[0]+"_vs_"+groups[1] + "-" + vt + ".nii.gz")

