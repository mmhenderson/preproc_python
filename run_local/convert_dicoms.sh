#!/bin/bash
#SBATCH --partition=tarrq
#SBATCH --gres=gpu:0
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --open-mode=append
#SBATCH --output=./sbatch_output/output-%A-%x-%u.out 
#SBATCH --time=8-00:00:00

# change this path
ROOT=/user_data/mmhender/

cd ${ROOT}preproc_code/PreprocScripts

inpath=${ROOT}data_UW/DataRaw/S01/Session1

echo ${inpath}

module load matlab-9.7

# echo ${PATH}
dcm2niix_path=${ROOT}preproc_code/dcm2niix
PATH=${PATH}:${dcm2niix_path}
echo ${PATH}
# echo ${dcm2niix_path}

echo starting matlab...
matlab -nodisplay -nodesktop -nosplash -r "convert_dicoms_forSiemens('"$inpath"'); exit"
echo finished


