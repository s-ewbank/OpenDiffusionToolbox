rootdir=/scratch/groups/rairan/NODDI/25-01-02/batch_output/

cd $rootdir
mkdir filt_results
subdirs=$(ls)

for d in $subdirs
do
cp $d/filt_log.csv filt_results/${d}_filt_log.csv
done