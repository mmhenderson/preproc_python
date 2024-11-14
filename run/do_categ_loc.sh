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
cd /lab_data/hendersonlab/preproc_code/

# subs=(S01 S02)
subs=(S02)

# module load fsl-6.0.3
module load freesurfer-7.1.0
source /user_data/mmhender/SetUpFreeSurfer.sh
# must have the above line for the subjects dir to be correct

# newer FSL version, seems to work better
# this is the same setup as my "newfsl" alias
FSLDIR=/lab_data/hawk/fsl-6.0.7.11
PATH=/lab_data/hawk/fsl-6.0.7.11/bin:$PATH
source /lab_data/hawk/fsl-6.0.7.11/etc/fslconf/fsl.sh

for sub in ${subs[@]}
do

    python3 -c 'from preproc_python import analyze_categ_loc; analyze_categ_loc.analyze_loc("'${sub}'")'  

done