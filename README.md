1. Prerequesite to using this toolkit, you should install ANTs, FSL, and 3 Singularity containers: Microstructure Diffusion Toolbox (MDT), MIRACL, and DSI-Studio. You can do this by running commands such as:
   ```
   singularity pull miracl_latest.sif library://aiconslab/miracl/miracl:latest
   singularity pull dsi_latest.sif library://dsistudio/dsistudio:latest
   ```
2. Download the toolkit.
3. Use the MAKE_CONFIG.sh script to make a config file by entering
   ```
   bash MAKE_CONFIG.sh
   ```
   and answering all subsequent prompts.
4. Use the script to fit DTI/NODDI metrics with
   ```
   bash RUN_SCRIPTS.sh fit
   ```
5. Then register with
   ```
   bash RUN_SCRIPTS.sh reg
   ```
6. Then quantify signal within atlas regions using:
   ```
   bash RUN_SCRIPTS.sh roi
   ```
7. Then quantify signal using voxel-based analysis with:
   ```
   bash RUN_SCRIPTS.sh vba
   ```
7. Finally, do tractography analysis with:
   ```
   bash RUN_SCRIPTS.sh tract
   ```
