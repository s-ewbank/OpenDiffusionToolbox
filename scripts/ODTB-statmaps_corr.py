###############################################
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
from scipy.stats import ttest_ind, pearsonr, spearmanr
from scipy.ndimage import label, generate_binary_structure
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
plt.rc('text', usetex=False)
rcParams['font.sans-serif'] = ['Arial']
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score, classification_report
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests
import argparse


###############################################
# INITIALIZE FILE, DIRECTORY NAMESff
###############################################
parser = argparse.ArgumentParser(description="Script for generating stat maps.")

parser.add_argument("--rootdir", required=True, help="root directory path")
parser.add_argument("--mask", required=True, help="path to mask")
parser.add_argument("--outdir", required=True, help="output directory path")
parser.add_argument("--avg_vol_dir", required=True, help="path to average volumes")
parser.add_argument("--data", required=True, help="path to scalar data file for correlation")
parser.add_argument("--organism", required=True, help="organism")
parser.add_argument("--vt", required=True, help="name of volume type (e.g. FA, MD)")
parser.add_argument("--design", required=True, help="statistical design; baseline, delta, or longitudinal")

args = parser.parse_args()
rootdir = args.rootdir
mask_path = args.mask
output_dir = args.outdir
avg_vol_dir = args.avg_vol_dir
corr_data_file = args.data
organism = args.organism
vt = args.vt
design = args.design

smooth=True

if organism=="mouse":
    cut_coords=np.arange(-3.5,3.5,1.5)
    fwhm=0.1
elif organism=="rat":
    cut_coords=np.arange(-13,6,2.5)
    fwhm=0.5
elif organism=="human":
    cut_coords=np.arange(-60,60,20)
    fwhm=5

corr_data_df=pd.read_csv(corr_data_file)
groups=sorted(corr_data_df["group"].unique())
if len(groups)>1:
    print(f"WARNING - multiple groups detected in correlation data file; all will be grouped together: {groups}")

def downsample_image(img, downsample_factor):
    original_affine = img.affine
    target_shape = tuple([int(d / downsample_factor) for d in img.shape])
    
    # Modify the diagonal elements of the affine matrix to scale by the downsample factor
    # This effectively increases the voxel size
    new_affine = original_affine.copy()
    new_affine[:3, :3] *= downsample_factor
    
    # Resample the image using the modified affine matrix and new target shape
    resampled_img = nilearn.image.resample_img(img, target_affine=new_affine, target_shape=target_shape, interpolation="nearest")
    
    return resampled_img

def threshold_and_cluster(pmap_img, tmap_img, voxel_threshold, voxel_cluster_count):
    voxel_size = pmap_img.header.get_zooms()[:3]
    voxel_volume = np.prod(voxel_size)
    min_region_size_mm3 = voxel_volume * voxel_cluster_count
    
    # Apply voxel threshold
    thresholded_map = nilearn.image.math_img(f"img < {voxel_threshold}", img=pmap_img)
    thresholded_data = thresholded_map.get_fdata()
    
    # Define connectivity structure
    s = generate_binary_structure(3, 2)
    
    # Label connected components
    labeled_array, n_labels = label(thresholded_data, structure=s)
    label_sizes = np.bincount(labeled_array.ravel())
    
    # Create a filter to keep labels with sizes above the threshold
    area_filter = label_sizes >= voxel_cluster_count
    area_filter[0] = False  # Ignore background label
    
    # Mask out smaller clusters
    filtered_data = area_filter[labeled_array]  # This creates a mask of 0s and 1s
    final_data = filtered_data * thresholded_data  # Apply the mask to original thresholded data
    filtered_labels = labeled_array * area_filter[labeled_array]
    
    # Create new Nifti image
    final_labels_img = nilearn.image.new_img_like(pmap_img, final_data)
    
    # Fetch the T-map data
    tmap_data = tmap_img.get_fdata()
    
    n_pos_clust = 0
    n_neg_clust = 0

    # Iterate over each cluster and check if its average T-value is positive or negative
    for i in range(1, len(area_filter)):  # starting from 1 to skip background
        if area_filter[i]:
            cluster_mask = filtered_labels == i
            t_values = tmap_data[cluster_mask]
            
            if np.mean(t_values) > 0:
                n_pos_clust += 1
            else:
                n_neg_clust += 1
                
    if n_labels == 0:
        print("No clusters identified!")
        final_map = np.zeros_like(tmap_data)
        final_map = nib.Nifti1Image(final_map, affine=pmap_img.affine, header=pmap_img.header)
    else:
        final_map = nilearn.image.math_img("tmap * labels", tmap=tmap_img, labels=final_labels_img)
    print(f"identified {n_pos_clust} positive clusters and {n_neg_clust} negative clusters")
    
    return final_map, final_labels_img, n_pos_clust, n_neg_clust




#groups_dict=load_groups(all_image_names)
#groups_dict=adjust_ntxgrps(groups_dict)



###############################################
# GETTING SUBJ IDS AND PLOTTING DIFFS AND T SCORES
###############################################
if design=="baseline":
    mask_data = np.squeeze(nib.load(mask_path).get_fdata())
    
    #cmap="PiYG_r"
    #cmap="BrBG_r"
    cmap="PuOr_r"
    
    voltypes=[vt]
    
    group=sorted(corr_data_df["group"].unique())
    
    mask_loaded=0
    
    clust_n = []
    
    for voltype in voltypes:
        imgs=[]
        imgs_data=[]
        subjs=0
        corr_data=[]
        img_names=list(corr_data_df["filename"])
        
        voltype_count=0
        
        print(f"Analyzing {voltype} for {group}")
    
        for subj in range(len(img_names)):
            subj_id=img_names[subj].split("_")[0]
            #try:
            subj_pre_dir=img_names[subj]
            subj_pre_vol=[rootdir+"/"+subj_pre_dir+"_output_reg_output/"+i for i in os.listdir(rootdir+"/"+subj_pre_dir+"_output_reg_output") if "_"+voltype+"_reg.nii.gz" in i][0]
            subjs+=1
            corr_data_i=corr_data_df[corr_data_df["filename"]==subj_pre_dir]["data"].to_numpy()[0]
            #except:
            #    print("Subj " + subj_id + " had some issue, so skipping.")
            #    continue
            if mask_loaded==0:
                mask=nilearn.image.new_img_like(nib.load(subj_pre_vol), mask_data)
                mask_loaded=1
                affine = nib.load(subj_pre_vol).affine
                header = nib.load(subj_pre_vol).header.copy()
                
            if smooth==True:
                subj_pre_vol=nilearn.image.smooth_img(subj_pre_vol, fwhm)
            else:
                subj_pre_vol=subj_pre_vol
                
            masked_img = nilearn.image.math_img("mask*pre",pre=subj_pre_vol,mask=mask)
            imgs.append(masked_img)
            imgs_data.append(masked_img.get_fdata())
            corr_data.append(corr_data_i)
            
        corrmap_name=f"corrmap_baseline_corr_{voltype}group-{group}"
        
        
        print("Computing correlation for volume " + voltype)
        del imgs
        del subj_pre_vol
        del masked_img
        imgs_data_nan = {}
        imgs_data_nan = np.where(imgs_data == 0, np.nan, imgs_data)
        #del imgs_data
        
        
        corrmap=np.zeros_like(imgs_data_nan[0])
        pmap=np.zeros_like(imgs_data_nan[0])
        n_slices=pmap.shape[0]
        y_dim=pmap.shape[1]
        z_dim=pmap.shape[2]
        
        for slice in range(n_slices):
            im_slice=[vol[slice,:,:] for vol in imgs_data_nan]
            corr_dat_slice=[np.ones_like(im_slice[0])*c for c in corr_data]
            result=pearsonr(im_slice,corr_dat_slice,axis=0)
            corrmap[slice,:,:]=result.statistic
            pmap[slice,:,:]=result.pvalue
        
        #for slice in range(n_slices):
        #    print(f"slice {slice+1} of {n_slices}")
        #    for y in range(y_dim):
        #        for z in range(z_dim):
        #            im_line=np.array([vol[slice,y,z] for vol in imgs_data])
        #            result=spearmanr(im_line,corr_data,nan_policy="omit")
        #            corrmap[slice,y,z]=result.statistic
        #            pmap[slice,y,z]=result.pvalue
        
        del result
        del imgs_data_nan
        
        #pmap_flat = pmap.ravel()
        #adjusted_p = multipletests(pmap_flat, alpha=0.1, method='fdr_bh')[1]
        #adjusted_pmap = adjusted_p.reshape(pmap.shape)
        
        corrmap_img = nib.Nifti1Image(corrmap, affine=affine, header=header)
        corrmap_img.to_filename(output_dir+"/"+corrmap_name+".nii.gz")
        pmap_img = nib.Nifti1Image(pmap, affine=affine, header=header)
        print(np.min(pmap_img))
        #adjusted_pmap_img = nib.Nifti1Image(adjusted_pmap, affine=affine, header=header)
        
        nilearn.plotting.plot_img(corrmap_img,title="Correlation Map Baseline {voltype} - {group}",cmap=cmap,
            draw_cross=False,annotate=False,vmin=-6,vmax=6,output_file=output_dir+f"/correlation_map_baseline_{voltype}-{group}.png",black_bg=True,
            display_mode="y",cut_coords=cut_coords,colorbar=True)
        
        print("Done computing t-test for volume " + voltype + ", doing clustering")
        
        voxel_threshold = 0.05
        cluster_size = 200
        final_tmap_img, labels_img, n_pos_clust, n_neg_clust = threshold_and_cluster(pmap_img, corrmap_img, voxel_threshold, cluster_size)
        
        clust_n.append(pd.DataFrame({"pos": [n_pos_clust], "neg": [n_neg_clust]},index=[str(f"clust_{voltype}group-{group}-corr")]))
        
        template_img = nib.load(avg_vol_dir+"/avg_"+voltype+"_masked.nii.gz")
        
        
        print("Saving result for " + voltype + " in output dir as " + corrmap_name+".nii.gz and "+corrmap_name+".png")
        
        final_tmap_img.to_filename(output_dir + "/thresholded_" + corrmap_name + ".nii.gz")
        
        # custom cmap
        print("plotting stat map")
        nilearn.plotting.plot_stat_map(final_tmap_img, bg_img=template_img, title=f"Thresh Corrmap Baseline {voltype} {group} ({subjs})", display_mode="y",annotate=False,
                                       threshold=voxel_threshold, cut_coords=cut_coords,cmap=cmap,black_bg=True,
                                       output_file=output_dir+"/thresholded_"+corrmap_name+".png",dim=1,vmax=1,vmin=-1)
        # If using pop, you can do: "/thresholded_"+corrmap_name+"_pop"+str(pop_number)+".png"
        
        del corrmap
        del pmap
        del pmap_img
        del corrmap_img
    
    final_clust_df = pd.concat(clust_n)
    final_clust_df.to_csv(output_dir+"/clusters_"+voltype+"_"+groups[0]+ "_v_"+groups[1]+".csv")
    