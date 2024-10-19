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
generate_tmaps=int(sys.argv[5])
do_classifier=int(sys.argv[6])
vt=sys.argv[7]
delta_or_baseline=int(sys.argv[8])
smooth=True

# Current bads
#Removing 10 images because they exceed 'bad frame' percent threshold of 5%:
#  - A3_Pre_-_240311183637_output_miracl_output
#  - B2_Post_-_240315100939_output_miracl_output
#  - M23_Ntx_Ketamine_Pre_-_231022143802_output_miracl_output
#  - M31_Ketamine_Post_-_231115120515_output_miracl_output
#  - M28_Ntx_Saline_Post_-_231025171718_output_miracl_output
#  - M17_Ntx_Ketamine_Post_-_231020123958_output_miracl_output
#  - M11_Saline_Post_-_230929121443_output_miracl_output
#  - M13_Ketamine_Pre_-_230928135856_output_miracl_output
#  - M18_Ntx_Saline_Pre_-_231019123325_output_miracl_output
#  - M8_Ketamine_Post_-_230921162519_output_miracl_output


groups_dict={}

groups_dict["F_sal"]=[['F1_Saline_Post_-_230511170358_output_miracl_output', 'F1_Saline_Pre_-_230510161927_output_miracl_output'],
['F10_Saline_Post_-_231116160102_output_miracl_output', 'F10_Saline_Pre_-_231115144642_output_miracl_output'],
['F12_Saline_Post_-_231201151809_output_miracl_output', 'F12_Saline_Pre_-_231130141524_output_miracl_output'],
['F14_Saline_Post_-_231206130048_output_miracl_output', 'F14_Saline_Pre_-_231205121554_output_miracl_output'],
['F16_Saline_Post_-_231206143436_output_miracl_output', 'F16_Saline_Pre_-_231205140531_output_miracl_output'],
['F18_Saline_Post_-_240612124123_output_miracl_output', 'F18_Saline_Pre_-_240611122412_output_miracl_output'],
['F3_Saline_Post_-_230511185433_output_miracl_output'],
['F5_Saline_Post_-_230516171413_output_miracl_output', 'F5_Saline_Pre_-_230515172819_output_miracl_output'],
['F6_Saline_Post_-_231031154459_output_miracl_output', 'F6_Saline_Pre_-_231030151331_output_miracl_output'],
['F8_Saline_Post_-_231115130042_output_miracl_output', 'F8_Saline_Pre_-_231114115523_output_miracl_output']]

groups_dict["M_sal"]=[['M1_Saline_Post_-_230503154113_output_miracl_output', 'M1_Saline_Pre_-_230502141952_output_miracl_output'],
['M10_Saline_Post_-_230929112838_output_miracl_output', 'M10_Saline_Pre_-_230928104652_output_miracl_output'],
#['M11_Saline_Post_-_230929121443_output_miracl_output', 'M11_Saline_Pre_-_230928115435_output_miracl_output'],
['M14_Saline_Post_-_230929144614_output_miracl_output', 'M14_Saline_Pre_-_230928144623_output_miracl_output'],
['M30_Saline_Post_-_231115111718_output_miracl_output', 'M30_Saline_Pre_-_231114101308_output_miracl_output'],
['M32_Saline_Post_-_240626133834_output_miracl_output', 'M32_Pre_-_240625130010_output_miracl_output'],
['M34_Saline_Post_-_240626151500_output_miracl_output', 'M34_Pre_-_240625144219_output_miracl_output'],
['M36_Saline_Post_-_240626164724_output_miracl_output', 'M36_Saline_Pre_-_240625162055_output_miracl_output'],
['M6_Saline_Post_-_230516162559_output_miracl_output', 'M6_Saline_Pre_-_230515153953_output_miracl_output'],
['M9_Saline_Post_-_230921175429_output_miracl_output', 'M9_Saline_Pre_-_230920173827_output_miracl_output']]


groups_dict["F_ket"]=[['F11_Ketamine_Post_-_231116164718_output_miracl_output', 'F11_Ketamine_Pre_-_231115153837_output_miracl_output'],
['F13_Ketamine_Post_-_231201160856_output_miracl_output', 'F13_Ketamine_Pre_-_231130150424_output_miracl_output'],
['F15_Ketamine_Post_-_231206134814_output_miracl_output', 'F15_Ketamine_Pre_-_231205132041_output_miracl_output'],
['F17_Ketamine_Post_-_231207181615_output_miracl_output', 'F17_Ketamine_Pre_-_231206152252_output_miracl_output'],
['F19_Ketamine_Post_-_240612133057_output_miracl_output', 'F19_Ketamine_Pre_-_240611131933_output_miracl_output'],
['F2_Ketamine_Post_-_230511180228_output_miracl_output', 'F2_Ketamine_Pre_-_230510171827_output_miracl_output'],
['F20_Ketamine_Post_-_240612141639_output_miracl_output', 'F20_Ketamine_Pre_-_240611140858_output_miracl_output'],
['F4_Ketamine_Post_-_230511194356_output_miracl_output'],
['F7_Ketamine_Post_-_231031163058_output_miracl_output', 'F7_Ketamine_Pre_-_231030160709_output_miracl_output'],
['F9_Ketamine_Post_-_231116151024_output_miracl_output', 'F9_Ketamine_Pre_-_231115135700_output_miracl_output']]

groups_dict["M_ket"]=[['M12_Ketamine_Post_-_230929130748_output_miracl_output', 'M12_Ketamine_Pre_-_230928124516_output_miracl_output'],
#['M13_Ketamine_Post_-_230929135734_output_miracl_output', 'M13_Ketamine_Pre_-_230928135856_output_miracl_output'],
['M15_Ketamine_Post_-_231004150847_output_miracl_output', 'M15_Ketamine_Pre_-_231003150421_output_miracl_output'],
['M3_Ketamine_Post_-_230503165331_output_miracl_output', 'M3_Ketamine_Pre_-_230502153909_output_miracl_output'],
#['M31_Ketamine_Post_-_231115120515_output_miracl_output', 'M31_Ketamine_Pre_-_231114110604_output_miracl_output'],
['M33_Ketamine_Post_-_240626142640_output_miracl_output', 'M33_Pre_-_240625135513_output_miracl_output'],
['M35_Ketamine_Post_-_240626160000_output_miracl_output', 'M35_Ketamine_Pre_-_240625153202_output_miracl_output'],
['M5_Ketamine_Post_-_230516153444_output_miracl_output', 'M5_Ketamine_Pre_-_230515144803_output_miracl_output'],
['M7_Ketamine_Pre_-_230515163507_output_miracl_output']]#,
#['M8_Ketamine_Post_-_230921162519_output_miracl_output', 'M8_Ketamine_Pre2_-_230920162549_output_miracl_output']]

groups_dict["M_ntxsal"]=[['M16_Ntx_Saline_Post_-_231020114653_output_miracl_output', 'M16_Ntx_Saline_Pre_-_231019104806_output_miracl_output'],
#['M18_Ntx_Saline_Post_-_231020132350_output_miracl_output', 'M18_Ntx_Saline_Pre_-_231019123325_output_miracl_output'],
['M20_Ntx_Saline_Post_-_231020151108_output_miracl_output', 'M20_Ntx_Saline_Pre_-_231019142309_output_miracl_output'],
['M22_Ntx_Saline_Post_-_231023135135_output_miracl_output', 'M22_Ntx_Saline_Pre_-_231022134714_output_miracl_output'],
['M24_Ntx_Saline_Post_-_231023152108_output_miracl_output', 'M24_Ntx_Saline_Pre_-_231022153049_output_miracl_output'],
['M26_Ntx_Saline_Post_-_231025154350_output_miracl_output'],
#['M28_Ntx_Saline_Post_-_231025171718_output_miracl_output','M28_Ntx_Saline_-_231024165055_output_miracl_output'],
['M37_Ntx_Saline_Post_-_240628130321_output_miracl_output', 'M37_Ntx_Saline_Pre_-_240627123821_output_miracl_output']]

groups_dict["M_ntxket"]=[#['M17_Ntx_Ketamine_Post_-_231020123958_output_miracl_output', 'M17_Ntx_Ketamine_Pre_-_231019113534_output_miracl_output'],
['M19_Ntx_Ketamine_Post_-_231020142121_output_miracl_output', 'M19_Ntx_Ketamine_Pre_-_231019132925_output_miracl_output'],
['M21_Ntx_Ketamine_Post_-_231023130225_output_miracl_output', 'M21_Ntx_Ketamine_Pre_-_231022130109_output_miracl_output'],
#['M23_Ntx_Ketamine_Post_-_231023143630_output_miracl_output', 'M23_Ntx_Ketamine_Pre_-_231022143802_output_miracl_output'],
['M25_Ntx_Ketamine_Post_-_231025144831_output_miracl_output','M25_Ntx_Ketamine_-_231024142414_output_miracl_output'],
['M27_Ntx_Ketamine_Post_-_231025163043_output_miracl_output','M27_Ntx_Ketamine_-_231024160047_output_miracl_output'],
['M29_Ntx_Ketamine_Post_-_231025180325_output_miracl_output','M29_Ntx_Ketamine_-_231024173840_output_miracl_output'],
['M38_Ntx_Ketamine_Post_-_240628135334_output_miracl_output', 'M38_Ntx_Ketamine_Pre_-_240627132738_output_miracl_output']]

groups_dict["SAPAP_ket"]=[["A1_Post_-_240314141001_output_miracl_output","A1_Pre_-_240311161724_output_miracl_output"],
["B1_Post_-_240315090933_output_miracl_output","B1_Pre_-_240312161728_output_miracl_output"],
#["B2_Post_-_240315100939_output_miracl_output","B2_Pre_-_240312171134_output_miracl_output"],
["C1_Post_-_240425101429_output_miracl_output","C1_Pre_-_240422124816_output_miracl_output"],
["C3_Post_-_240425120702_output_miracl_output","C3_Pre_-_240422143307_output_miracl_output"],
["D3_Post_-_240627160403_output_miracl_output","D3_Pre_-_240624143605_output_miracl_output"]]

groups_dict["SAPAP_sal"]=[["A2_Post_-_240314151334_output_miracl_output","A2_Pre_-_240311174322_output_miracl_output"],
#["A3_Post_-_240314160732_output_miracl_output","A3_Pre_-_240311183637_output_miracl_output"],
["B3_Post_-_240315105400_output_miracl_output","B3_Pre_-_240312181033_output_miracl_output"],
["C2_Post_-_240425110942_output_miracl_output","C2_Pre_-_240422133749_output_miracl_output"],
["D1_Post_-_240627141802_output_miracl_output","D1_Pre_-_240624124601_output_miracl_output"],
["D2_Post_-_240627151159_output_miracl_output","D2_Pre_-_240624134154_output_miracl_output"]]

groups_dict["WT_ket"]=[["E3_Post_-_240711173647_output_miracl_output","E3_Pre_-_240708174437_output_miracl_output"],
["E1_Post_-_240711155906_output_miracl_output","E1_Pre_-_240708160722_output_miracl_output"],
["E2_Post_-_240711164620_output_miracl_output","E2_Pre_-_240708165714_output_miracl_output"]]

groups_dict["WT_sal"]=[["G3_Post_-_240726152215_output_miracl_output","G3_Pre_-_240723162349_output_miracl_output"],
["G1_Post_-_240726134657_output_miracl_output","G1_Pre_-_240723143453_output_miracl_output"],
["G2_Post_-_240726143259_output_miracl_output","G2_Pre_-_240723152530_output_miracl_output"]]

baseline_groups_dict={}

M_keys=["M_ket","M_sal"]#,"M_ntxsal","M_ntxket"]
M=[]
for key in M_keys:
    M=M+[i[1] for i in groups_dict[key] if len(i)>1]
baseline_groups_dict["M"]=M

F_keys=["F_sal","F_ket"]
F=[]
for key in F_keys:
    F=F+[i[1] for i in groups_dict[key] if len(i)>1]
baseline_groups_dict["F"]=F

SAPAP_keys=["SAPAP_ket","SAPAP_sal"]
SAPAP=[]
for key in SAPAP_keys:
    SAPAP=SAPAP+[i[1] for i in groups_dict[key] if len(i)>1]
baseline_groups_dict["SAPAP"]=SAPAP

WT_keys=["WT_ket","WT_sal"]
WT=[]
for key in WT_keys:
    WT=WT+[i[1] for i in groups_dict[key] if len(i)>1]
baseline_groups_dict["WT"]=WT

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
        final_map = []
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
        
        groups=["M_ntxket","M_ntxsal"]
        
        voltypes=[vt]
        
        mask_loaded=0
        
        clust_n = []
        
        for voltype in voltypes:
            diff_imgs_data={}
            for group in groups:
                diff_imgs=[]
                image_pairs=groups_dict[group]
                
                voltype_count=0
                
                diff_imgs_data[group]=[]
                print("Analyzing " + voltype + " for " + group)
            
                for subj in range(len(image_pairs)):
                    subj_id=image_pairs[subj][0].split("_")[0]
                    if len(image_pairs[subj])==2:
                        try:
                            subj_pre_dir=image_pairs[subj][1]
                            subj_pre_vol=[rootdir+"/"+subj_pre_dir+"/"+i for i in os.listdir(rootdir+"/"+subj_pre_dir) if "_"+voltype+"_reg.nii.gz" in i][0]
                            subj_post_dir=image_pairs[subj][0]
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
                        fwhm=0.1
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
                    nilearn.plotting.plot_img(masked_diff_img,title="Delta "+voltype+" " +subj_id+" - "+group,cmap="seismic",
                        draw_cross=False,annotate=False,vmin=vmin,vmax=vmax,output_file=output_dir+"/individual_diffs/delta_"+voltype+"_"+subj_id+"_"+group+".png",
                        display_mode="y",cut_coords=np.arange(-13,-1,1.5))
                
                #Compute mean diff
                mean_diff=nilearn.image.mean_img(diff_imgs)
                
                #Save diff image to nifti
                mean_diff.to_filename(output_dir+"/delta_"+voltype+"_"+group+".nii.gz")
                
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
                nilearn.plotting.plot_img(mean_diff,title="Delta "+voltype+" " +group,cmap="seismic",
                    draw_cross=False,annotate=False,vmin=vmin,vmax=vmax,output_file=output_dir+"/delta_"+voltype+"_"+group+".png",
                    display_mode="y",cut_coords=np.arange(-13,-1,1.5),colorbar=True)
                    
            tmap_name="tmap_delta_"+voltype+"group1-"+groups[0]+"group2-"+groups[1]
            
            
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
                result=ttest_ind(group0_slice,group1_slice, nan_policy='omit')
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
            #adjusted_pmap_img = nib.Nifti1Image(adjusted_pmap, affine=affine, header=header)
            
            nilearn.plotting.plot_img(tmap_img,title="T-Score Delta "+voltype+" " +groups[0]+ " vs. "+groups[1],cmap="seismic",
                draw_cross=False,annotate=False,vmin=-6,vmax=6,output_file=output_dir+"/"+tmap_name+".png",black_bg=True,
                display_mode="y",cut_coords=np.arange(-13,-1,1.5),colorbar=True)
            
            print("Done computing t-test for volume " + voltype + ", doing clustering")
            
            voxel_threshold = 0.05
            cluster_size = 200
            final_tmap_img, labels_img, n_pos_clust, n_neg_clust = threshold_and_cluster(pmap_img, tmap_img, voxel_threshold, cluster_size)
            
            clust_n.append(pd.DataFrame({"pos": [n_pos_clust], "neg": [n_neg_clust]},index=[str("clust_"+voltype+"group1-"+groups[0]+"group2-"+groups[1])]))
            
            template_img = nib.load(avg_vol_dir+"/avg_"+voltype+"_masked.nii.gz")
            
            
            print("Saving result for " + voltype + " in output dir as " + tmap_name+".nii.gz and "+tmap_name+".png")
            
            final_tmap_img.to_filename(output_dir + "/thresholded_" + tmap_name + ".nii.gz")
            
            # custom cmap
            nilearn.plotting.plot_stat_map(final_tmap_img, bg_img=template_img, title="Thresh T-Score Delta "+voltype+" " +groups[0]+ " vs. "+groups[1], display_mode="y",annotate=False,
                                           threshold=voxel_threshold, cut_coords=np.arange(-13,-1,1.5),cmap="seismic",black_bg=True,
                                           output_file=output_dir+"/thresholded_"+tmap_name+".png",dim=1,vmax=12)
            # If using pop, you can do: "/thresholded_"+tmap_name+"_pop"+str(pop_number)+".png"
            
            del tmap
            del pmap
            del pmap_img
            del tmap_img
        
        final_clust_df = pd.concat(clust_n)
        final_clust_df.to_csv(output_dir+"/clusters_"+voltype+"_"+groups[0]+ "_v_"+groups[1]+".csv")
        
elif delta_or_baseline==1:
    if generate_tmaps==1:
        mask_data = np.squeeze(nib.load(mask_path).get_fdata())
        
        groups=["F","M"]
        cmap="PiYG_r"
        
        #groups=["SAPAP","WT"]
        #cmap="BrBG_r"
        
        voltypes=[vt]
        
        mask_loaded=0
        
        clust_n = []
        
        
        for voltype in voltypes:
            masked_img_data={}
            output_fnames={}
            for group in groups:
                output_fnames[group]=[]
                print("Analyzing " + voltype + " for " + group)
                for subj in range(len(baseline_groups_dict[group])):
                    subj_id=baseline_groups_dict[group][subj].split("_")[0]
                    subj_dir=baseline_groups_dict[group][subj]
                    subj_vol=[rootdir+"/"+subj_dir+"/"+i for i in os.listdir(rootdir+"/"+subj_dir) if "_"+voltype+"_reg.nii.gz" in i][0]
                    if mask_loaded==0:
                        mask=nilearn.image.new_img_like(nib.load(subj_vol), mask_data)
                        mask_loaded=1
                    if smooth==True:
                        fwhm=0.1
                        subj_vol=nilearn.image.smooth_img(subj_vol, fwhm)
                    else:
                        subj_vol=nib.load(subj_vol)
                    masked_img = nilearn.image.math_img("mask*vol",vol=subj_vol,mask=mask)
                    masked_img_data = masked_img.get_fdata()
                    masked_img_data_nan = np.where(masked_img_data == 0, np.nan, masked_img_data)
                    output_fname=str(output_dir+"/tmp/"+group+"_"+subj_id)
                    output_fnames[group].append(output_fname)
                    np.save(output_fname, masked_img_data_nan)
                    del masked_img
                    del masked_img_data
                    del masked_img_data_nan
                    del subj_vol
                    print("Done with subj " + subj_id)

                    
        tmap_name="tmap_"+voltype+"group1-"+groups[0]+"group2-"+groups[1]
        tmap=np.zeros_like(mask_data)
        pmap=np.zeros_like(mask_data)
        n_slices=pmap.shape[0]
        
        for slice in range(n_slices):
            masked_img_data_nan={}
            for group in groups:
                masked_img_data_nan[group]=[]
                for subj in output_fnames[group]:
                    subj_data=np.load(subj+".npy")
                    masked_img_data_nan[group].append(subj_data[slice,:,:])
                    del subj_data
            print("Making t/p maps for slice " + str(slice) + " of " + str(n_slices))
            result=ttest_ind(masked_img_data_nan[groups[0]],masked_img_data_nan[groups[1]], nan_policy='omit')
            tmap[slice,:,:]=result.statistic
            pmap[slice,:,:]=result.pvalue
            
            del result
            del masked_img_data_nan
            
            print("Done computing t-test for volume " + voltype + ", slice " + str(slice) + " of " + str(n_slices))
            #del diff_imgs_data
            
            #pmap_flat = pmap.ravel()
            #adjusted_p = multipletests(pmap_flat, alpha=0.1, method='fdr_bh')[1]
            #adjusted_pmap = adjusted_p.reshape(pmap.shape)
            
        affine = mask.affine
        header = mask.header.copy()
        tmap_img = nib.Nifti1Image(tmap, affine=affine, header=header)
        tmap_img.to_filename(output_dir+"/"+tmap_name+".nii.gz")
        pmap_img = nib.Nifti1Image(pmap, affine=affine, header=header)
        #adjusted_pmap_img = nib.Nifti1Image(adjusted_pmap, affine=affine, header=header)
        
        nilearn.plotting.plot_img(tmap_img,title="T-Score Delta "+voltype+" " +groups[0]+ " vs. "+groups[1],cmap=cmap,
            draw_cross=False,annotate=False,vmin=-6,vmax=6,output_file=output_dir+"/"+tmap_name+".png",black_bg=True,
            display_mode="y",cut_coords=np.arange(-13,-1,1.5),colorbar=True)
        
        print("Done computing t-test for volume " + voltype + ", doing clustering")
        
        voxel_threshold = 0.05
        cluster_size = 200
        final_tmap_img, labels_img, n_pos_clust, n_neg_clust = threshold_and_cluster(pmap_img, tmap_img, voxel_threshold, cluster_size)
        
        clust_n.append(pd.DataFrame({"pos": [n_pos_clust], "neg": [n_neg_clust]},index=[str("clust_"+voltype+"group1-"+groups[0]+"group2-"+groups[1])]))
        
        template_img = nib.load(avg_vol_dir+"/avg_"+voltype+"_masked.nii.gz")
        
        
        print("Saving result for " + voltype + " in output dir as " + tmap_name+".nii.gz and "+tmap_name+".png")
        
        final_tmap_img.to_filename(output_dir + "/baseline_thresholded_" + tmap_name + ".nii.gz")
        
        # custom cmap
        nilearn.plotting.plot_stat_map(final_tmap_img, bg_img=template_img, title="Thresh T-Score Baseline "+voltype+" " +groups[0]+ " vs. "+groups[1], display_mode="y",annotate=False,
                                       threshold=voxel_threshold, cut_coords=np.arange(-13,-1,1.5),cmap=cmap,black_bg=True,
                                       output_file=output_dir+"/baseline_thresholded_"+tmap_name+".png",dim=1,vmax=6)
        # If using pop, you can do: "/thresholded_"+tmap_name+"_pop"+str(pop_number)+".png"
        
        del tmap
        del pmap
        del pmap_img
        del tmap_img
        
        final_clust_df = pd.concat(clust_n)
        final_clust_df.to_csv(output_dir+"/baseline_clusters_"+voltype+"_"+groups[0]+ "_v_"+groups[1]+".csv")