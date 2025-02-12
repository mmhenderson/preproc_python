#!/bin/bash

# making an output filename that has my timestamp
outstr=$(basename $0)_$(date '+%Y-%m-%d_%H:%M:%S')
{
# source /home/lab/miniconda3/envs/fmrienv/bin/activate
# source activate fmrienv
echo $CONDA_PREFIX
# conda list


# where your copy of preprocessing code lives
code_path=/home/lab/hendersonlab/code_featsynth/preproc_code/
cd $code_path

# can change this path - where your data lives
data_path=/home/lab/hendersonlab/data_featsynth

sub=S03
sub_FS=AD

# python3 -c 'from preproc_python import make_roi_volumes; make_roi_volumes.proc_retino_annot('"${sub_FS}"')'
# python3 -c 'from preproc_python import make_roi_volumes; make_roi_volumes.proc_floc_annot('"${sub_FS}"')'
echo python3 -c 'from preproc_python import make_roi_volumes; make_roi_volumes.make_vols_func("'${sub}'","'${sub_FS}'")'

python3 -c 'from preproc_python import make_roi_volumes; make_roi_volumes.make_vols_func("'${sub}'","'${sub_FS}'")'

} 2>&1 | tee -a ./sbatch_output/output_${outstr}.txt >/dev/null
# IF you want output to terminal too, take out the /dev/null part
