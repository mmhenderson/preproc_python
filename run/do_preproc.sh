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


sub=S04
sub_retino=AE

debug=0

# # module load fsl-6.0.3
# module load freesurfer-7.1.0
# source /user_data/mmhender/SetUpFreeSurfer.sh
# # must have the above line for the subjects dir to be correct

# # newer FSL version, seems to work better
# # this is the same setup as my "newfsl" alias
# FSLDIR=/lab_data/hawk/fsl-6.0.7.11
# PATH=/lab_data/hawk/fsl-6.0.7.11/bin:$PATH
# source /lab_data/hawk/fsl-6.0.7.11/etc/fslconf/fsl.sh


python3 -c 'from preproc_python import preprocessing; preprocessing.run_preproc("'${sub}'","'${sub_retino}'","'${debug}'")'  


} 2>&1 | tee -a ./sbatch_output/output_${outstr}.txt >/dev/null
# IF you want output to terminal too, take out the /dev/null part

