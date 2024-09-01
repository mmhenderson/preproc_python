#!/bin/bash
#SBATCH --partition=tarrq
#SBATCH --gres=gpu:0
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --open-mode=append
#SBATCH --output=./sbatch_output/output-%A-%x-%u.out 
#SBATCH --time=8-00:00:00

source /user_data/mmhender/fmriproc_env/bin/activate

# change this path
ROOT=/user_data/mmhender/

cd ${ROOT}preproc_code/

sub=S01
sess_list=(Sess2)

# sub=S02
# sess_list=(Sess1)

# echo $PATH

module load fsl-6.0.3

for sess in ${sess_list[@]}
do

    inpath=${ROOT}data_featsynth/DataRaw/${sub}/${sess}
    echo ${inpath}

    python3 -c 'from preproc_python import unwarping; unwarping.do_topup("'${inpath}'")'  

done