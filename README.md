# OpenDiffusionToolbox
## Description
Toolkit for implementing dMRI model fitting (DTI, NODDI), registration, and ROI- and voxel-based analysis in mice, rats, or humans using an HPC like Stanford's Sherlock.
## Installation
1. Download this repository and the atlases that go with it. The atlases can be installed from the following link:
   ```
    wget --no-check-certificate -O /opt/atlases.zip -d "https://www.dropbox.com/scl/fo/1iy3woqtvl7f5vhet33jb/AMwaSBtm2M2O8mvfDIcwXn0?rlkey=wqege46bidvnjty6pz0x1w8nk&st=b4elsg99&dl=1"
   ```
2. I have installed *almost all the dependencies for this toolbox in a Singularity container for which the Singularity.def makefile is provided. The Singularity container can be made as follows, and should be usable without any additional exporting of paths or mounting:
   ```
   singularity build ODTB.sif /path/to/OpenDiffusionToolbox/singularity_recipe/Singularity11.def
   ```
*If you wish to use MDT (e.g., for NODDI model fitting), you should separately use the container recipe provided by [akhanf](https://github.com/akhanf/mdt-singularity) to make another container to run in conjunction with this, since it has the appropriate setup for OpenCL kernel management.
## Usage
### Before you start
Before beginning, you should format your data such that data for each scan/observation is in its own folder containing the following: (1) a 4D raw diffusion MRI file in NIFTI format; (2) a bvecs file with suffix .bvec, and (3) a bvals file with suffix .bval. For example:
   ```
   subject1/
   ├── 4d_dmri_img.nii.gz     # 4D raw dMRI image (.nii.gz)
   ├── bvals.bval             # bval file (.bval)
   ├── bvecs.bvec             # bvec file (.bvec)
   
   subject2/
   ├── 4d_dmri_img.nii.gz
   ├── bvals.bval
   ├── bvecs.bvec
   
   subject3/
   ├── 4d_dmri_img.nii.gz
   ├── bvals.bval
   ├── bvecs.bvec
   
   ```
### Processing
Once you have this, you may proceed to the following steps. Making the config and doing the "fit" and "reg" steps are prerequesite to all subsequent steps in the "Analysis" section (i.e., roi, tract, vba_compare)
1. First, use the MAKE_CONFIG.sh script to make a config file with the following (and note that the parent directory of these subdirectories is what you will enter as root directory in the config file):
   ```
   bash MAKE_CONFIG.sh
   ```
   and answering all subsequent prompts.
2. Use the script to fit DTI/NODDI metrics with
   ```
   bash RUN_SCRIPTS.sh --config <path/to/config> --step fit
   ```
3. Then register with
   ```
   bash RUN_SCRIPTS.sh --config <path/to/config> --step reg
   ```
### Analysis
After doing the "fit" and "reg" steps, you may move on to doing analysis! Note that the analyses are based on the groups CSV file specified in the config as "groups_file", and a dummy version of this file needing only unique timepoint and condition labels is made during creation of the config file.
4. Quantify signal within atlas regions using:
   ```
   bash RUN_SCRIPTS.sh --config <path/to/config> --step roi
   ```
5. Get full-study average volumes with:
   ```
   bash RUN_SCRIPTS.sh --config <path/to/config> --step vba_avg
   ```
6. Do tractography analysis with:
   ```
   bash RUN_SCRIPTS.sh --config <path/to/config> --step tract
   ```
7. Do tbss with:
   ```
   bash RUN_SCRIPTS.sh --config <path/to/config> --step tbss
   ```
8. Finally, run voxel-based group comparisons using vba_compare with the added --design flag set to a valid option (baseline, delta, longitudinal for full volume or baseline_tbss, delta_tbss, longitudinal_tbss for TBSS voxels). The comparison will use the group assignments designated in the groups CSV file indicated in the config file. An example is as follows:
   ```
   bash RUN_SCRIPTS.sh --config <path/to/config> --step vba_compare --design delta
   ```