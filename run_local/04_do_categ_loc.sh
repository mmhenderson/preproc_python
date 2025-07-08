#!/bin/bash

# making an output filename that has my timestamp
outstr=$(basename $0)_$(date '+%Y-%m-%d_%H:%M:%S')
{
# source /home/lab/miniconda3/envs/fmrienv/bin/activate
# source activate fmrienv
echo $CONDA_PREFIX

# where your copy of preprocessing code lives
code_path=/home/lab/hendersonlab/code_featsynth/preproc_code/
cd $code_path

# subs=(S01 S02)
subs=(S03)
subs_FS=(AD)

# module load fsl-6.0.3
# module load freesurfer-7.1.0
# source /user_data/mmhender/SetUpFreeSurfer.sh
# # must have the above line for the subjects dir to be correct

# # newer FSL version, seems to work better
# # this is the same setup as my "newfsl" alias
# FSLDIR=/lab_data/hawk/fsl-6.0.7.11
# PATH=/lab_data/hawk/fsl-6.0.7.11/bin:$PATH
# source /lab_data/hawk/fsl-6.0.7.11/etc/fslconf/fsl.sh

for si in "${!subs[@]}"
# for sub in ${subs[@]}
do

    sub=${subs[$si]}
    sub_fs=${subs_FS[$si]}
    
    python3 -c 'from preproc_python import analyze_categ_loc; analyze_categ_loc.analyze_loc("'${sub}'", "'${sub_fs}'")'  

done


} 2>&1 | tee -a ./sbatch_output/output_${outstr}.txt >/dev/null
# IF you want output to terminal too, take out the /dev/null part

