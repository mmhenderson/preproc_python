#!/bin/bash
#SBATCH --partition=tarrq
#SBATCH --gres=gpu:0
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --open-mode=append
#SBATCH --output=./sbatch_output/output-%A-%x-%u.out 
#SBATCH --time=8-00:00:00

source /user_data/mmhender/fmriproc_env/bin/activate


# source activate fmrienv
echo $CONDA_PREFIX

# where your copy of preprocessing code lives
code_path=/lab_data/hendersonlab/code_featsynth/preproc_code/
cd $code_path

# subj=(11)
subj=(8)

for ss in ${subj[@]}
do

    python3 -c 'from preproc_python import process_behav; process_behav.proc("'${ss}'")'  

done


# } 2>&1 | tee -a ./sbatch_output/output_${outstr}.txt >/dev/null
# # IF you want output to terminal too, take out the /dev/null part

