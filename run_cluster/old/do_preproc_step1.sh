#!/bin/bash
#SBATCH --partition=tarrq
#SBATCH --gres=gpu:0
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --open-mode=append
#SBATCH --output=./sbatch_output/output-%A-%x-%u.out 
#SBATCH --time=8-00:00:00

echo $SLURM_JOBID
echo $SLURM_NODELIST


source /user_data/mmhender/fmriproc_env/bin/activate

# change this path
cd /user_data/mmhender/preproc_code/

sub=S01B
sub_retino=AA
debug=0

module load fsl-6.0.3
module load freesurfer-7.1.0
source /user_data/mmhender/SetUpFreeSurfer.sh
# must have the above line for the subjects dir to be correct

python3 -c 'from code import preprocessing; preprocessing.run_preproc_step1("'${sub}'","'${sub_retino}'","'${debug}'")'  
