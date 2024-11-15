#!/bin/bash
#SBATCH -p normal
#SBATCH -c 1
#SBATCH --cpus-per-task=4
#SBATCH --time 8:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=snewbank@stanford.edu
#SBATCH --output=/home/groups/rairan/NODDI/slurm_outputs//slurm-%j.out
#

rootdir="/scratch/groups/rairan/NODDI/24-08-31/batch_volumes_raw_240831/"
preregdir="/scratch/groups/rairan/NODDI/24-08-31/batch_volumes_raw_240831/batch_output/"

cd $rootdir
subdirs=$(ls)
n_subdirs=`find . -mindepth 1 -maxdepth 1 -type d | wc -l`
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Identified $n_subdirs subdirectories in root directory for DSI fiber tracking in $rootdir"
echo " "
count=1

for d in $subdirs
do
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running DSI tractography script on subdir $count called $d."
    bash /home/groups/rairan/NODDI/experimental/sne-dsi.sh $rootdir $preregdir $d
    count=$(($count+1))
done