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
1. First, use the MAKE_CONFIG.sh script to make a config file by entering
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
4. Then quantify signal within atlas regions using:
   ```
   bash RUN_SCRIPTS.sh --config <path/to/config> --step roi
   ```
5. Then quantify signal using voxel-based analysis with:
   ```
   bash RUN_SCRIPTS.sh --config <path/to/config> --step vba
   ```
6. Finally, do tractography analysis with:
   ```
   bash RUN_SCRIPTS.sh --config <path/to/config> --step tract
   ```