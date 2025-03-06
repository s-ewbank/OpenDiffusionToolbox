import nibabel as nib
import sys
import numpy as np

qa_in=sys.argv[1]
mask=sys.argv[2]
qa_out=sys.argv[3]

qa = nib.load(qa_in)
qa_data = qa.get_fdata()

mask = nib.load(mask)

qa_data_flipped = np.flip(qa_data, axis=0)
qa_align = nib.Nifti1Image(qa_data_flipped, mask.affine, mask.header)
print("Saving qa align file at "+qa_out)
nib.save(qa_align, qa_out)