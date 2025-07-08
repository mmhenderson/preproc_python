#!/bin/bash

# making an output filename that has my timestamp
outstr=$(basename $0)_$(date '+%Y-%m-%d_%H-%M-%S')
{
# source /home/lab/miniconda3/envs/fmrienv/bin/activate
source activate /home/lab/miniconda3/envs/fmrienv

# where your copy of preprocessing code lives
code_path=/home/lab/hendersonlab/code_featsynth/preproc_code/
cd $code_path

sub=S06
sub_retino=AG

debug=0

python3 -c 'from preproc_python import copy_ret_files; copy_ret_files.copy_all("'${sub}'","'${sub_retino}'","'${debug}'")'  


} 2>&1 | tee -a ./sbatch_output/output_${outstr}.txt >/dev/null
# IF you want output to terminal too, take out the /dev/null part
