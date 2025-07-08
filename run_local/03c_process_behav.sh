##!/bin/bash

# making an output filename that has my timestamp
outstr=$(basename $0)_$(date '+%Y-%m-%d_%H:%M:%S')
{
# source /home/lab/miniconda3/envs/fmrienv/bin/activate
# source activate fmrienv
echo $CONDA_PREFIX

# where your copy of preprocessing code lives
code_path=/home/lab/hendersonlab/code_featsynth/preproc_code/
cd $code_path

subj=(4)

for ss in ${subj[@]}
do

    python3 -c 'from preproc_python import process_behav; process_behav.proc("'${ss}'")'  

done


} 2>&1 | tee -a ./sbatch_output/output_${outstr}.txt >/dev/null
# IF you want output to terminal too, take out the /dev/null part

