#!/bin/bash
#SBATCH -c 4            # Number of CPUS requested. If omitted, the default is 1 CPU.
#SBATCH --mem=10240     # Memory requested in megabytes. If omitted, the default is 1024 MB.
#SBATCH -t 0-20:0:0      # How long will your job run for? If omitted, the default is 3 hours.

matlab -nodesktop -nosplash -r 'run(/global/home/hpc3586/JE_packages/pca_je/pca_fmri.m)'