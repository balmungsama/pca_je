#!/bin/bash
#SBATCH -c 4            # Number of CPUS requested. If omitted, the default is 1 CPU.
#SBATCH --mem=10240     # Memory requested in megabytes. If omitted, the default is 1024 MB.
#SBATCH -t 0-20:0:0      # How long will your job run for? If omitted, the default is 3 hours.

## grab user input 

while getopts t:o:p:b:f: option; do
	case "${option}" in
		t) top_dir=${OPTARG}
			 top_dir="top_dir='$top_dir'";;
		o) output=${OPTARG}
			 output="output='$output'";;
		p) pipe=${OPTARG}
			 pipe="pipe=$pipe";;
		b) nboot=${OPTARG}
			 nboot="nboot=$nboot";;
		f) filters=${OPTARG}
			 filters="filters=$filters"
	esac
done

matlab_cmd="run('/global/home/hpc3586/JE_packages/pca_je/pca_on_fmri.m')"
matlab_cmd="$top_dir;$output;$pipe;$nboot;$filters;$matlab_cmd"

matlab -nodesktop -nosplash -r $matlab_cmd