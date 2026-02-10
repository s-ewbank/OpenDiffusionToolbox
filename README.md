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
2. Use the script to fit DTI/NODDI metrics by setting the step flag to "fit":
   ```
   bash RUN_SCRIPTS.sh --config <path/to/config> --step fit
   ```
3. Then run a two-step (linear and nonlinear) registration to an atlas with step flag set to "reg":
   ```
   bash RUN_SCRIPTS.sh --config <path/to/config> --step reg
   ```
### ROI-based analysis
After doing the "fit" and "reg" steps, you may move on to doing analysis! As a first step, we can extract the signal in all atlas regions.
4. Quantify signal within atlas regions with step flag set to "roi":
   ```
   bash RUN_SCRIPTS.sh --config <path/to/config> --step roi
   ```
### Spatial voxel-wise analysis
A next step is to do spatial voxel-wise analysis, either through TBSS or through whole-brain voxel-wise analysis. Note that the analyses are based on the groups CSV file specified in the config as "groups_file", and a dummy version of this file needing only unique timepoint and condition labels is made during creation of the config file. For executing voxel-wise statistical comparisons, I have commands to use FSL randomise as well as my own custom Python-based voxel-wise, cluster size thresholded comparisons.
5. Get full-study average volumes with step flag set to "vba_avg" (you will need this for voxel-wise steps!):
   ```
   bash RUN_SCRIPTS.sh --config <path/to/config> --step vba_avg
   ```
6a. Extract tbss skeletons for all volumes with step flag set to "tbss_1":
   ```
   bash RUN_SCRIPTS.sh --config <path/to/config> --step tbss
   ```
6b. Then, execute tract-based spatial statistical comparison with FSL randomise (for group files as defined in groups .csv file specified in config) by setting step flag to "tbss_2":
   ```
   bash RUN_SCRIPTS.sh --config <path/to/config> --step tbss
   ```
7a. For whole-brain voxel-wise analysis with FSL randomise for statistical computations, set step flag to "vba" and add design flag set to "baseline" (for 2 group, 1 timepoint comparisons), "delta" (for 2 group, 2 timepoint comparisons), or "correlate" (for 1 group with a variable to correlate added to the groups file in a new column):
   ```
   bash RUN_SCRIPTS.sh --config <path/to/config> --step vba_randomise --design <baseline/delta/correlate>
   ```
7b. Alternatively, run whole-brain voxel-wise analysis with my custom cluster-based, voxel-wise analysis approach with step flag set to "vba_compare" and design set to "baseline", "delta", or "longitudinal".
   ```
   bash RUN_SCRIPTS.sh --config <path/to/config> --step vba_compare --design <baseline/delta/longitudinal>
   ```
### Tractography-based analysis
A next step is to do tractography based analysis. These implementations (particularly probabalistic tractography) have not been as thoroughly tested as the other method implementations (e.g., registration, ROI-based, voxel-wise)
8. Do DSI-based deterministic tractography analysis with step flag set to "tract":
   ```
   bash RUN_SCRIPTS.sh --config <path/to/config> --step tract
   ```
9. Do FSL-based probabalistic tractography analysis with step flag set to "probtrackx":
   ```
   bash RUN_SCRIPTS.sh --config <path/to/config> --step tract
   ```