#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --gres=gpu:0
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --open-mode=append
#SBATCH --output=./sbatch_output/output-%A-%x-%u.out 
#SBATCH --time=8-00:00:00

source /user_data/mmhender/fmriproc_env/bin/activate

echo $CONDA_PREFIX
# conda list


# where your copy of preprocessing code lives
code_path=/lab_data/hendersonlab/code_featsynth/preproc_code/
cd $code_path

# can change this path - where your data lives
data_path=/lab_data/hendersonlab/data_featsynth

sub=S08
# sess_list=(Sess1 Sess2 Sess3 Sess4)
sess_list=(Sess2)

for sess in ${sess_list[@]}
do

    inpath=${data_path}/DataRaw/${sub}/${sess}
    echo ${inpath}

    # python3 -c 'from preproc_python import prep_data; prep_data.process_tar_raw("'${inpath}'")'  
    # python3 -c 'from preproc_python import prep_data; prep_data.organize_functional("'${inpath}'")'  
    python3 -c 'from preproc_python import prep_data; prep_data.get_run_info("'${inpath}'")'  

done

# } 2>&1 | tee -a ./sbatch_output/output_${outstr}.txt >/dev/null
# # IF you want output to terminal too, take out the /dev/null part
