# DiffusionMRIToolkit
## Description
Toolkit for implementing dMRI model fitting (DTI, NODDI), registration, and ROI- and voxel-based analysis in mice, rats, or humans using an HPC like Stanford's Sherlock.
## Installation
1. Prerequesite to using this toolkit, you should install ANTs, FSL, and 3 Singularity containers: Microstructure Diffusion Toolbox (MDT), MIRACL, and DSI-Studio. The MDT container should be installed as instructed by [akhanf](https://github.com/akhanf/mdt-singularity), and the MIRACL and DSI-Studio containers can be built as follows:
   ```
   singularity pull miracl_latest.sif library://aiconslab/miracl/miracl:latest
   singularity pull dsi_latest.sif library://dsistudio/dsistudio:latest
   ```
2. Download the toolkit.

## Usage
1. First, use the MAKE_CONFIG.sh script to make a config file by entering
   ```
   bash MAKE_CONFIG.sh
   ```
   and answering all subsequent prompts.
2. Use the script to fit DTI/NODDI metrics with
   ```
   bash RUN_SCRIPTS.sh fit
   ```
3. Then register with
   ```
   bash RUN_SCRIPTS.sh reg
   ```
4. Then quantify signal within atlas regions using:
   ```
   bash RUN_SCRIPTS.sh roi
   ```
5. Then quantify signal using voxel-based analysis with:
   ```
   bash RUN_SCRIPTS.sh vba
   ```
6. Finally, do tractography analysis with:
   ```
   bash RUN_SCRIPTS.sh tract
   ```
