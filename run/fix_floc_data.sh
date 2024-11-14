#!/bin/bash
#SBATCH --partition=tarrq
#SBATCH --gres=gpu:0
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --open-mode=append
#SBATCH --output=./sbatch_output/output-%A-%x-%u.out 
#SBATCH --time=8-00:00:00

# this is where a copy of floc code is located
ROOT=/user_data/mmhender/
cd ${ROOT}preproc_code/fLoc_mod

# this part gets changed for your experiment
data_path=/user_data/mmhender/data_UW/DataBehavior/S01/Session1/

module load matlab-9.7

echo ${data_path}

# NOTE this won't work unless PTB is installed. so really you should run it on the scanner 
# laptop, or any personal machine where PTB is installed.
echo starting matlab...
matlab -nodisplay -nodesktop -nosplash -r "convert_floc_data('"$data_path"'); exit"
echo finished


