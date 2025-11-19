#!/bin/bash
#
### !!! CHANGE !!! the email address to your drexel email
#SBATCH --mail-user=ar953@drexel.edu
### !!! CHANGE !!! the account - you need to consult with the professor
#SBATCH --account=eces450650prj

### select number of nodes (usually you need only 1 node)
#SBATCH --nodes=1
### select number of tasks per node
#SBATCH --ntasks=1
### select number of cpus per task (you need to tweak this when you run a multi-thread program)
#SBATCH --cpus-per-task=1
### request 1 hour of wall clock time (if you request less time, you can wait for less time to get your job run by the system, you need to have a good esitmation of the run time though).
#SBATCH --time=06:00:00
### memory size required per node (this is important, you also need to estimate a upper bound)
#SBATCH --mem=16GB
### select the partition "edu" (this is the default partition but you can change according to your application)
#SBATCH --partition=edu


### Whatever modules you used (e.g. picotte-openmpi/gcc)
### must be loaded to run your code.
### Add them below this line.

# Define variables
MOTUS_DOCKER_IMAGE=docker://quay.io/biocontainers/motus:3.0.1--pyhdfd78af_0
LOCAL_SIF_NAME=motus_3.0.1.sif


# containerdir=/motus_tutorial
echo "Pulling Singularity image ${LOCAL_SIF_NAME}..."
singularity pull --name motus_3.0.1.sif docker://quay.io/biocontainers/motus:3.0.1--pyhdfd78af_0
singularity exec --bind ~/motus_db:/usr/local/share/motus-3.0.1/db_mOTU ${LOCAL_SIF_NAME} bash motus_profiler_t7.sh


### The commands you want to run in this job script (run a python script, run a certain software with inputs and outpus etc.).
echo "Executing mOTUs analysis script inside container..."
singularity exec ${LOCAL_SIF_NAME} bash motus_profiler_t7.sh

echo "job finished"
