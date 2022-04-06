#!/bin/bash
#SBATCH --output /scratch/lab_boztug/akamnev/1/JobLogs/%j.log   # log file location (both stdout and stderr), %j stands for unique job ID
#SBATCH --job-name=CPProfilingAPID
#SBATCH --partition=shortq                                       # job queue where the job is submitted to
#SBATCH --ntasks=1                                                # 1 task
# Optional parameters
#SBATCH --cpus-per-task=1                                         # 1 task on 1 CPU
#SBATCH --mem-per-cpu=6000                                        # 1 task on 1 CPU 1*4Gb of memory
#SBATCH --time=2:59:00                                           # Job time is max 30 minutes
#SBATCH --error /scratch/lab_boztug/akamnev/1/JobLogs/%j.err  # error log file location (stderr), %j stands for unique job ID
# SBATCH --nodes=1                                                 # number of nodes
#SBATCH --exclude=i001

#add path to modules for calling?
source /etc/profile.d/modules.sh #fixes bag with sbatch --export skipping some jobs

# Get modules needed to run CellProfiler
module load gcc/4.8.2 python/2.7.12 java/jdk/1.8.0/102 cellprofiler/3.1.8

date

echo Batch $iter

# Define variables
batchFileName=Batch_data.h5
batchFilePath=/scratch/lab_boztug/akamnev/Plate/$batchFileName
firstImage=$((9*$iter - 8))
lastImage=$(($firstImage + 8))
echo Processing images $firstImage to $lastImage

# Run CellProfiler in batch mode
cellprofiler -p $batchFilePath -c -r -f $firstImage -l $lastImage -o /scratch//lab_boztug/akamnev/Plate/Results/batch_$iter/

date
