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
from scipy.stats import ttest_ind
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

###############################################
# INITIALIZE FILE, DIRECTORY NAMES
###############################################
rootdir=sys.argv[1]
mask_path=sys.argv[2]
output_dir=sys.argv[3]
avg_vol_dir=sys.argv[4]
groups_file=sys.argv[5]
organism=sys.argv[6]
generate_tmaps=int(sys.argv[7])
do_classifier=int(sys.argv[8])
vt=sys.argv[9]
delta_or_baseline=int(sys.argv[10])
smooth=True

if organism=="mouse":
    cut_coords=np.arange(-9,-2,1.5)
    fwhm=0.1
elif organism=="rat":
    cut_coords=np.arange(-13,6,2.5)
    fwhm=0.5
elif organism=="human":
    cut_coords=np.arange(-60,60,20)
    fwhm=5

groups_df=pd.read_csv(groups_file)
groups=groups_df["group"].unique()

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
if delta_or_baseline==0:
    if generate_tmaps==1:
        mask_data = np.squeeze(nib.load(mask_path).get_fdata())
        
        groups=sorted(list(timepts_dict.keys()))
        
        groupnames={}
        
        for group in groups:
            print("Group "+str(group)+" contains:")
            [print(i) for i in timepts_dict[group]]
            groupnames[group]=input("What should this group be called? ")
        
        voltypes=[vt]
        
        mask_loaded=0
        
        clust_n = []
        
        for voltype in voltypes:
            diff_imgs_data={}
            for group in groups:
                diff_imgs=[]
                image_pairs=timepts_dict[group]
                
                voltype_count=0
                
                diff_imgs_data[group]=[]
                print("Analyzing " + voltype + " for " + str(group))
            
                for subj in range(len(image_pairs)):
                    subj_id=image_pairs[subj][0].split("_")[0]
                    if len(image_pairs[subj])==2:
                        try:
                            subj_pre_dir=image_pairs[subj][0]
                            subj_pre_vol=[rootdir+"/"+subj_pre_dir+"/"+i for i in os.listdir(rootdir+"/"+subj_pre_dir) if "_"+voltype+"_reg.nii.gz" in i][0]
                            subj_post_dir=image_pairs[subj][1]
                            subj_post_vol=[rootdir+"/"+subj_post_dir+"/"+i for i in os.listdir(rootdir+"/"+subj_post_dir) if "_"+voltype+"_reg.nii.gz" in i][0]
                        except:
                            print("Subj " + subj_id + " had some issue, so skipping.")
                            continue
                    else:
                        print("Subj " + subj_id + " did not have a pre and or post image, so skipping.")
                        continue
                    if mask_loaded==0:
                        mask=nilearn.image.new_img_like(nib.load(subj_pre_vol), mask_data)
                        mask_loaded=1
                        
                    if smooth==True:
                        subj_post_vol=nilearn.image.smooth_img(subj_post_vol, fwhm)
                        subj_pre_vol=nilearn.image.smooth_img(subj_pre_vol, fwhm)
                    else:
                        subj_post_vol=subj_post_vol
                        subj_pre_vol=subj_pre_vol
                    masked_diff_img = nilearn.image.math_img("mask*(post-pre)",post=subj_post_vol,pre=subj_pre_vol,mask=mask)
                    diff_imgs.append(masked_diff_img)
                    
                    diff_imgs_data[group].append(masked_diff_img.get_fdata())
                    
                    
                    print("Done with subj " + subj_id)
                    
                    #Make a png
                    if voltype=="FA":
                        vmin=-0.25
                        vmax=0.25
                    elif voltype=="NDI":
                        vmin=-0.4
                        vmax=0.4
                    elif voltype=="ODI":
                        vmin=-0.4
                        vmax=0.4
                    elif voltype=="MD":
                        vmin=-0.0004
                        vmax=0.0004
                    nilearn.plotting.plot_img(masked_diff_img,title="Delta "+voltype+" " +subj_id+" - "+groupnames[group],cmap="seismic",
                        draw_cross=False,annotate=False,vmin=vmin,vmax=vmax,output_file=output_dir+"/individual_diffs/delta_"+voltype+"_"+subj_id+"_"+groupnames[group]+".png",
                        display_mode="y",cut_coords=cut_coords)
                
                #Compute mean diff
                mean_diff=nilearn.image.mean_img(diff_imgs)
                
                #Save diff image to nifti
                mean_diff.to_filename(output_dir+"/delta_"+voltype+"_"+groupnames[group]+".nii.gz")
                
                #Make a png
                if voltype=="FA":
                    vmin=-0.25
                    vmax=0.25
                elif voltype=="NDI":
                    vmin=-0.4
                    vmax=0.4
                elif voltype=="ODI":
                    vmin=-0.4
                    vmax=0.4
                elif voltype=="MD":
                    vmin=-0.0004
                    vmax=0.0004
                nilearn.plotting.plot_img(mean_diff,title="Delta "+voltype+" " +groupnames[group],cmap="seismic",
                    draw_cross=False,annotate=False,vmin=vmin,vmax=vmax,output_file=output_dir+"/delta_"+voltype+"_"+groupnames[group]+".png",
                    display_mode="y",cut_coords=cut_coords,colorbar=True)
                    
            tmap_name="tmap_delta_"+voltype+"group1-"+groupnames[groups[0]]+"group2-"+groupnames[groups[1]]
            
            
            print("Computing t-test for volume " + voltype)
            del diff_imgs
            del subj_post_vol
            del subj_pre_vol
            del masked_diff_img
            diff_imgs_data_nan = {}
            for group in groups:
                diff_imgs_data_nan[group] = np.where(diff_imgs_data[group] == 0, np.nan, diff_imgs_data[group])
            del diff_imgs_data
            
            
            tmap=np.zeros_like(diff_imgs_data_nan[groups[0]][0])
            pmap=np.zeros_like(diff_imgs_data_nan[groups[0]][0])
            n_slices=pmap.shape[0]
            
            for slice in range(n_slices):
                group0_slice=[vol[slice,:,:] for vol in diff_imgs_data_nan[groups[0]]]
                group1_slice=[vol[slice,:,:] for vol in diff_imgs_data_nan[groups[1]]]
                result=ttest_ind(group1_slice,group0_slice, nan_policy='omit')
                tmap[slice,:,:]=result.statistic
                pmap[slice,:,:]=result.pvalue
            
            del result
            del diff_imgs_data_nan
            
            print("Done computing t-test for volume " + voltype + ", saving individual group tmaps")
            #del diff_imgs_data
            
            #pmap_flat = pmap.ravel()
            #adjusted_p = multipletests(pmap_flat, alpha=0.1, method='fdr_bh')[1]
            #adjusted_pmap = adjusted_p.reshape(pmap.shape)
            
            affine = mean_diff.affine
            header = mean_diff.header.copy()
            tmap_img = nib.Nifti1Image(tmap, affine=affine, header=header)
            tmap_img.to_filename(output_dir+"/"+tmap_name+".nii.gz")
            pmap_img = nib.Nifti1Image(pmap, affine=affine, header=header)
            print(np.min(pmap_img))
            #adjusted_pmap_img = nib.Nifti1Image(adjusted_pmap, affine=affine, header=header)
            
            nilearn.plotting.plot_img(tmap_img,title="T-Score Delta "+voltype+" " +groupnames[groups[0]]+ " vs. "+groupnames[groups[1]],cmap="seismic",
                draw_cross=False,annotate=False,vmin=-6,vmax=6,output_file=output_dir+"/"+tmap_name+".png",black_bg=True,
                display_mode="y",cut_coords=cut_coords,colorbar=True)
            
            print("Done computing t-test for volume " + voltype + ", doing clustering")
            
            voxel_threshold = 0.05
            cluster_size = 200
            final_tmap_img, labels_img, n_pos_clust, n_neg_clust = threshold_and_cluster(pmap_img, tmap_img, voxel_threshold, cluster_size)
            
            clust_n.append(pd.DataFrame({"pos": [n_pos_clust], "neg": [n_neg_clust]},index=[str("clust_"+voltype+"group1-"+groupnames[groups[0]]+"group2-"+groupnames[groups[1]])]))
            
            template_img = nib.load(avg_vol_dir+"/avg_"+voltype+"_masked.nii.gz")
            
            
            print("Saving result for " + voltype + " in output dir as " + tmap_name+".nii.gz and "+tmap_name+".png")
            
            final_tmap_img.to_filename(output_dir + "/thresholded_" + tmap_name + ".nii.gz")
            
            # custom cmap
            print("plotting stat map")
            nilearn.plotting.plot_stat_map(final_tmap_img, bg_img=template_img, title="Thresh T-Score Delta "+voltype+" " +groupnames[groups[0]]+ " vs. "+groupnames[groups[1]], display_mode="y",annotate=False,
                                           threshold=voxel_threshold, cut_coords=cut_coords,cmap="seismic",black_bg=True,
                                           output_file=output_dir+"/thresholded_"+tmap_name+".png",dim=1,vmax=12)
            # If using pop, you can do: "/thresholded_"+tmap_name+"_pop"+str(pop_number)+".png"
            
            del tmap
            del pmap
            del pmap_img
            del tmap_img
        
        final_clust_df = pd.concat(clust_n)
        final_clust_df.to_csv(output_dir+"/clusters_"+voltype+"_"+groupnames[groups[0]]+ "_v_"+groupnames[groups[1]]+".csv")

if delta_or_baseline==1:
    if generate_tmaps==1:
        mask_data = np.squeeze(nib.load(mask_path).get_fdata())
        
        groups=sorted(list(timepts_dict.keys()))
        
        #cmap="PiYG_r"
        cmap="BrBG_r"
        
        groupnames={}
        
        for group in groups:
            print("Group "+str(group)+" contains:")
            [print(i) for i in timepts_dict[group]]
            groupnames[group]=input("What should this group be called? ")
        
        voltypes=[vt]
        
        mask_loaded=0
        
        clust_n = []
        
        for voltype in voltypes:
            imgs_data={}
            for group in groups:
                imgs=[]
                image_pairs=timepts_dict[group]
                
                voltype_count=0
                
                imgs_data[group]=[]
                print("Analyzing " + voltype + " for " + str(group))
            
                for subj in range(len(image_pairs)):
                    subj_id=image_pairs[subj][0].split("_")[0]
                    try:
                        subj_pre_dir=image_pairs[subj][0]
                        subj_pre_vol=[rootdir+"/"+subj_pre_dir+"/"+i for i in os.listdir(rootdir+"/"+subj_pre_dir) if "_"+voltype+"_reg.nii.gz" in i][0]
                    except:
                        print("Subj " + subj_id + " had some issue, so skipping.")
                        continue
                    if mask_loaded==0:
                        mask=nilearn.image.new_img_like(nib.load(subj_pre_vol), mask_data)
                        mask_loaded=1
                        
                    if smooth==True:
                        subj_pre_vol=nilearn.image.smooth_img(subj_pre_vol, fwhm)
                    else:
                        subj_pre_vol=subj_pre_vol
                    masked_img = nilearn.image.math_img("mask*pre",pre=subj_pre_vol,mask=mask)
                    imgs.append(masked_img)
                    imgs_data[group].append(masked_img.get_fdata())
                    
                    print("Done with subj " + subj_id)
                    
                    #Make a png
                    if voltype=="FA":
                        vmin=-0.25
                        vmax=0.25
                    elif voltype=="NDI":
                        vmin=-0.4
                        vmax=0.4
                    elif voltype=="ODI":
                        vmin=-0.4
                        vmax=0.4
                    elif voltype=="MD":
                        vmin=-0.0004
                        vmax=0.0004
                    #nilearn.plotting.plot_img(masked_img,title="Baseline "+voltype+" " +subj_id+" - "+groupnames[group],cmap="gray",black_bg=True,
                    #    draw_cross=False,annotate=False,vmin=vmin,vmax=vmax,output_file=output_dir+"/individual_diffs/baseline_"+voltype+"_"+subj_id+"_"+groupnames[group]+".png",
                    #    display_mode="y",cut_coords=cut_coords)
                
                #Compute mean diff
                mean_im=nilearn.image.mean_img(imgs)
                
                #Save diff image to nifti
                mean_im.to_filename(output_dir+"/baseline_"+voltype+"_"+groupnames[group]+".nii.gz")
                
                #Make a png
                if voltype=="FA":
                    vmin=-0.25
                    vmax=0.25
                elif voltype=="NDI":
                    vmin=-0.4
                    vmax=0.4
                elif voltype=="ODI":
                    vmin=-0.4
                    vmax=0.4
                elif voltype=="MD":
                    vmin=-0.0004
                    vmax=0.0004
                nilearn.plotting.plot_img(mean_im,title="Baseline "+voltype+" " +groupnames[group],cmap="gray",black_bg=True,
                    draw_cross=False,annotate=False,vmin=vmin,vmax=vmax,output_file=output_dir+"/baseline_"+voltype+"_"+groupnames[group]+".png",
                    display_mode="y",cut_coords=cut_coords,colorbar=True)
                    
            tmap_name="tmap_baseline_"+voltype+"group1-"+groupnames[groups[0]]+"group2-"+groupnames[groups[1]]
            
            
            print("Computing t-test for volume " + voltype)
            del imgs
            del subj_pre_vol
            del masked_img
            imgs_data_nan = {}
            for group in groups:
                imgs_data_nan[group] = np.where(imgs_data[group] == 0, np.nan, imgs_data[group])
            del imgs_data
            
            
            tmap=np.zeros_like(imgs_data_nan[groups[0]][0])
            pmap=np.zeros_like(imgs_data_nan[groups[0]][0])
            n_slices=pmap.shape[0]
            
            for slice in range(n_slices):
                group0_slice=[vol[slice,:,:] for vol in imgs_data_nan[groups[0]]]
                group1_slice=[vol[slice,:,:] for vol in imgs_data_nan[groups[1]]]
                result=ttest_ind(group1_slice,group0_slice, nan_policy='omit')
                tmap[slice,:,:]=result.statistic
                pmap[slice,:,:]=result.pvalue
            
            del result
            del imgs_data_nan
            
            print("Done computing t-test for volume " + voltype + ", saving individual group tmaps")
            
            #pmap_flat = pmap.ravel()
            #adjusted_p = multipletests(pmap_flat, alpha=0.1, method='fdr_bh')[1]
            #adjusted_pmap = adjusted_p.reshape(pmap.shape)
            
            affine = mean_im.affine
            header = mean_im.header.copy()
            tmap_img = nib.Nifti1Image(tmap, affine=affine, header=header)
            tmap_img.to_filename(output_dir+"/"+tmap_name+".nii.gz")
            pmap_img = nib.Nifti1Image(pmap, affine=affine, header=header)
            print(np.min(pmap_img))
            #adjusted_pmap_img = nib.Nifti1Image(adjusted_pmap, affine=affine, header=header)
            
            nilearn.plotting.plot_img(tmap_img,title="T-Score Baseline "+voltype+" " +groupnames[groups[0]]+ " vs. "+groupnames[groups[1]],cmap=cmap,
                draw_cross=False,annotate=False,vmin=-6,vmax=6,output_file=output_dir+"/"+tmap_name+".png",black_bg=True,
                display_mode="y",cut_coords=cut_coords,colorbar=True)
            
            print("Done computing t-test for volume " + voltype + ", doing clustering")
            
            voxel_threshold = 0.05
            cluster_size = 200
            final_tmap_img, labels_img, n_pos_clust, n_neg_clust = threshold_and_cluster(pmap_img, tmap_img, voxel_threshold, cluster_size)
            
            clust_n.append(pd.DataFrame({"pos": [n_pos_clust], "neg": [n_neg_clust]},index=[str("clust_"+voltype+"group1-"+groupnames[groups[0]]+"group2-"+groupnames[groups[1]])]))
            
            template_img = nib.load(avg_vol_dir+"/avg_"+voltype+"_masked.nii.gz")
            
            
            print("Saving result for " + voltype + " in output dir as " + tmap_name+".nii.gz and "+tmap_name+".png")
            
            final_tmap_img.to_filename(output_dir + "/thresholded_" + tmap_name + ".nii.gz")
            
            # custom cmap
            print("plotting stat map")
            nilearn.plotting.plot_stat_map(final_tmap_img, bg_img=template_img, title="Thresh T-Score Baseline "+voltype+" " +groupnames[groups[0]]+ " vs. "+groupnames[groups[1]], display_mode="y",annotate=False,
                                           threshold=voxel_threshold, cut_coords=cut_coords,cmap=cmap,black_bg=True,
                                           output_file=output_dir+"/thresholded_"+tmap_name+".png",dim=1,vmax=6)
            # If using pop, you can do: "/thresholded_"+tmap_name+"_pop"+str(pop_number)+".png"
            
            del tmap
            del pmap
            del pmap_img
            del tmap_img
        
        final_clust_df = pd.concat(clust_n)
        final_clust_df.to_csv(output_dir+"/clusters_"+voltype+"_"+groupnames[groups[0]]+ "_v_"+groupnames[groups[1]]+".csv")
        
