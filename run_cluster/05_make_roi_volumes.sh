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


# source activate fmrienv
echo $CONDA_PREFIX
# conda list


# where your copy of preprocessing code lives
code_path=/lab_data/hendersonlab/code_featsynth/preproc_code/
cd $code_path

# sub=S03
# sub_FS=AD

sub=S11
sub_FS=AL



# module load fsl-6.0.3
module load freesurfer-7.1.0
source /user_data/mmhender/SetUpFreeSurfer.sh
# must have the above line for the subjects dir to be correct

# newer FSL version, seems to work better
# this is the same setup as my "newfsl" alias
FSLDIR=/lab_data/hawk/fsl-6.0.7.11
PATH=/lab_data/hawk/fsl-6.0.7.11/bin:$PATH
source /lab_data/hawk/fsl-6.0.7.11/etc/fslconf/fsl.sh


# python3 -c 'from preproc_python import make_roi_volumes; make_roi_volumes.proc_retino_annot('"${sub_FS}"')'
# python3 -c 'from preproc_python import make_roi_volumes; make_roi_volumes.proc_floc_annot('"${sub_FS}"')'

# echo python3 -c 'from preproc_python import make_roi_volumes; make_roi_volumes.make_vols_func("'${sub}'","'${sub_FS}'")'

python3 -c 'from preproc_python import make_roi_volumes; make_roi_volumes.make_vols_func("'${sub}'","'${sub_FS}'")'
python3 -c 'from preproc_python import make_roi_volumes; make_roi_volumes.threshold_ROI_masks("'${sub}'")'


# } 2>&1 | tee -a ./sbatch_output/output_${outstr}.txt >/dev/null
# # IF you want output to terminal too, take out the /dev/null part
