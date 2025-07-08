#!/bin/bash

PATH="/home/lab/miniconda3/bin/"

echo $PATH
# source /home/lab/miniconda3/envs/fmrienv/bin/activate
source activate fmrienv
echo $CONDA_PREFIX
# conda list

python3 -c 'import os; print(os.environ.get("CONDA_DEFAULT_ENV")); print(os.environ.get("CONDA_PREFIX"))'

python3 -c 'help("modules"); exit'

python3 -c 'import numpy as np; print(np.random.normal(0,1,10)); exit'
