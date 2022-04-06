#.sh file to call iteratively cluster CP jobs:

batchFileName=Cluster_RunSingleJob.sh
for i in {1..312}; do
  sbatch --export=iter=$i --output /scratch//lab_boztug/akamnev/Plate/JobLogs/iter$i"_%j.log" --job-name=$i"_CPProfilingAPID" --error /scratch//lab_boztug/akamnev/Plate/JobLogs/iter$i"_%j.err" $batchFileName;
  sleep 3;
done
