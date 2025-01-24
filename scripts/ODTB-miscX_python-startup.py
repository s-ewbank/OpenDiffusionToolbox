import nibabel as nib

# Monkey patch nibabel
def patch_nibabel():
    def patched_get_data(self):
        return self.get_fdata()
    nib.Nifti1Image.get_data = patched_get_data
    print("nibabel.get_data() has been patched to get_fdata() globally.")

patch_nibabel()