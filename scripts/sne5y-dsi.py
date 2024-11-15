import nibabel as nib
import sys
import numpy as np

FA=sys.argv[1]
labels_in=sys.argv[2]
out_path=sys.argv[3]


FA=nib.load(FA)
aff=FA.affine
hd=FA.header
labels=nib.load(labels_in)
labels_data=labels.get_fdata()
labels_data=labels_data.astype(int)

abbrev_labels_data=np.zeros_like(labels_data)

rois={"Frontal_ctx": [1,3,4,5,6,33,25],
    "Temporal_lobe": [8,9,10,11,12,13,14,15,16],
    "Anterior_cingulate_area": [29,30],
    "Somatosensory_areas": [17,18],
    "Somatomotor_areas": [7],
    "Thalamus": [104,115],
    "Striatum": [116,117,105,106]}

#rois={"mPFC": [44, 972],
#    "Anterior_cingulate_area": [39, 48],
#    "Somatosensory_areas": [329, 337, 345, 353, 361, 369, 378],
#    "Somatomotor_areas": [993, 500, 500],
#    "Visual_areas": [402, 394, 409, 425, 385, 533],
#    "Thalamus": [239, 958, 1008, 1014, 51, 138, 444, 571, 406, 856, 864, 637, 709],
#    "Hippocampus": [375, 726, 10704, 909, 1080, 822],
#    "Retrosplenial_area": [254, 254, 879, 894, 886],
#    "Pallidum": [904, 809, 818, 826, 835],
#    "Striatum": [275, 242, 485, 493],
#    "Amygdala": [278]}
    
for r, roi in enumerate(list(rois.keys())):
    val=r+1
    with open(out_path+'/lut_abbrev.txt', 'a') as file:
        file.write(roi+", "+str(val)+"\n")
    for id in rois[roi]:
        mask = labels_data == id
        abbrev_labels_data = np.where(mask, val, abbrev_labels_data)
    

new_img=nib.Nifti1Image(labels_data,affine=aff,header=hd)
nib.save(new_img,out_path+"/atlas_warp_flirt_RPI_hd.nii.gz")

new_img2=nib.Nifti1Image(abbrev_labels_data,affine=aff,header=hd)
nib.save(new_img2,out_path+"/atlas_warp_abbrev.nii.gz")