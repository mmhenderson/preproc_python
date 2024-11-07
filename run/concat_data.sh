#!/bin/bash
#SBATCH --partition=tarrq
#SBATCH --gres=gpu:0
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --open-mode=append
#SBATCH --output=./sbatch_output/output-%A-%x-%u.out 
#SBATCH --time=8-00:00:00

source /user_data/mmhender/fmriproc_env/bin/activate

# where your copy of preprocessing code lives
code_path=/lab_data/hendersonlab/preproc_code/
cd $code_path

subj=(2)

for ss in ${subj[@]}
do

    python3 -c 'from preproc_python import concat_data; concat_data.get_data_concat_maintask("'${ss}'")'  

done
