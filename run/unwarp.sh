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

# can change this path - where your data lives
data_path=/home/lab/hendersonlab/data_featsynth

sub=S03
sess_list=(Sess1 Sess2 Sess3 Sess4)


# stuff commented out below was to make it work on mind cluster. should not
# need this on lab workstation, because bashrc is already set up right

# we need all the below stuff to make sure afni / fsl paths are right
# if errors come up in do_topup, its usually to do with these paths.

# # module load fsl-6.0.3
# module load freesurfer-7.1.0
# source /user_data/mmhender/SetUpFreeSurfer.sh
# # must have the above line for the subjects dir to be correct

# # newer FSL version, seems to work better
# # this is the same setup as my "newfsl" alias
# FSLDIR=/lab_data/hawk/fsl-6.0.7.11
# PATH=/lab_data/hawk/fsl-6.0.7.11/bin:$PATH
# source /lab_data/hawk/fsl-6.0.7.11/etc/fslconf/fsl.sh

# PATH=/lab_data/hawk/afni_21.1.10:$PATH
# echo $PATH


for sess in ${sess_list[@]}
do

    inpath=${data_path}/DataRaw/${sub}/${sess}
    echo ${inpath}

    python3 -c 'from preproc_python import unwarping; unwarping.do_topup("'${inpath}'")'  

done

} 2>&1 | tee -a ./sbatch_output/output_${outstr}.txt >/dev/null
# IF you want output to terminal too, take out the /dev/null part
