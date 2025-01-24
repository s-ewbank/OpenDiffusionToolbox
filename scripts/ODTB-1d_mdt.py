import nibabel as nib
import mdt
import sys

dwi=sys.argv[1]
mask=sys.argv[2]
prtcl=sys.argv[3]

# Monkey patch to redirect get_data() to get_fdata()
def patch_nibabel():
    def patched_get_data(self):
        return self.get_fdata()

    nib.Nifti1Image.get_data = patched_get_data
    print("nibabel.get_data() has been patched to get_fdata() globally.")

# Apply the patch
patch_nibabel()

input_data = mdt.load_input_data(dwi,prtcl,mask)
mdt.fit_model('NODDI', input_data, 'output/noddi_output_dir/')