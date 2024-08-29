import numpy as np
import os, sys
import subprocess
import pandas as pd
from datetime import datetime

# where my topup script is. 
# Other users can make a copy of the script and edit this path.
topup_script_path = '/user_data/mmhender/preproc_code/preproc_python/unwarp/run_topup_forSiemens'


def do_topup(inpath, do_moco = False):

    # Do unwarping (aka topup, aka distortion correction, aka topup correction...).
    # This step takes your raw functional data (nifti files), and uses your 
    # field maps to correct for inhomogeneities in the magnetic field.
    
    # calls "unwarp/run_topup_forSiemens"
    
    original_dir = os.getcwd()
    
    if do_moco:
        outstr = '_topup_moco'
        moco_str = ' -domoco'
    else:
        # usually leave this off, we'll do motion correction later.
        outstr = '_topup'
        moco_str = ''

    # This is a pd dataframe that we already made.
    # It has info about which runs we did in which order.
    # Includes both field maps and functional scans.
    info_fn = os.path.join(inpath, 'Niftis', 'run_info.csv')
    print('loading info from %s'%info_fn)
    sys.stdout.flush()
    run_info = pd.read_csv(info_fn)
   
    # find different kinds of runs 
    func_inds = [i for i in np.arange(run_info.shape[0]) if 'fmap' not in np.array(run_info['run_type'])[i]]
    
    # these are the field maps we will use
    fmap_fwd_inds = [i for i in np.arange(run_info.shape[0]) if 'fmap_PA' in np.array(run_info['run_type'])[i]]
    fmap_rev_inds = [i for i in np.arange(run_info.shape[0]) if 'fmap_AP' in np.array(run_info['run_type'])[i]]

    # looping over scans.
    for ii, fi in enumerate(func_inds):
        
        # full name of the file we're working with here.
        file_to_unwarp = np.array(run_info['nifti_fn'])[fi].split('.nii.gz')[0]
        
        # going into the folder holding my file of interest. 
        # all outputs get created here.
        go_to_path = os.path.dirname(file_to_unwarp)
        os.chdir(go_to_path)
        
        # this is the name of the file we're going to create. 
        # could call it something else if you want.
        # it gets spit out inside of whatever folder we run topup command from
        output_stem = file_to_unwarp.split('/')[-1] + outstr

        # checking if it's already done...if so then we should skip.
        fn_topup = os.path.join(go_to_path, output_stem + '.nii.gz')
        print(fn_topup)
        sys.stdout.flush()
        if os.path.exists(fn_topup):
            print('done with this file, skipping')
            sys.stdout.flush()
            continue
        
        # check if tmp folder is already here, and if so get rid of it.
        # this would happen if the function got stuck somewhere in the middle.
        if os.path.exists(os.path.join(go_to_path,'tmp')):
            print('removing %s...'%os.path.join(go_to_path,'tmp'))
            err = subprocess.call('rm -r %s'%os.path.join(go_to_path,'tmp'), shell=True)
             
        
        # decide which fmaps we should use for unwarping here.
        # if there is just one pair for the session, we'll always use that.
        # otherwise choosing closest in time.
        t_here = np.array(run_info['time_raw'])[fi]

        t_fmaps = [np.array(run_info['time_raw'])[x] for x in fmap_fwd_inds]
        closer_fmap_ind = np.argmin(np.abs(np.array(t_fmaps)-t_here))
        fmap_fwd_ind_use = fmap_fwd_inds[closer_fmap_ind]

        t_fmaps = [np.array(run_info['time_raw'])[x] for x in fmap_rev_inds]
        closer_fmap_ind = np.argmin(np.abs(np.array(t_fmaps)-t_here))
        fmap_rev_ind_use = fmap_rev_inds[closer_fmap_ind]

        # get actual nifti filenames for the field maps.
        # removing the .nii.gz extension here.
        fmap_fwd_fn = np.array(run_info['nifti_fn'])[fmap_fwd_ind_use].split('.nii.gz')[0]
        fmap_rev_fn = np.array(run_info['nifti_fn'])[fmap_rev_ind_use].split('.nii.gz')[0]

        print(fmap_fwd_fn, fmap_rev_fn, file_to_unwarp)
        sys.stdout.flush()
        
        

        # this is where we actually call the topup script.
        # making a linux command here which we will execute.
        cmd = topup_script_path + \
                ' -d1 ' + fmap_fwd_fn + \
                ' -d2 ' + fmap_rev_fn + \
                ' -i ' + file_to_unwarp + \
                ' -o ' + output_stem + \
                moco_str

        print(cmd)
        sys.stdout.flush()
        
        err = subprocess.call(cmd, shell=True)

        if err:
            print('Error during topup function for file %d, %s'%(fi, file_to_unwarp))
            sys.stdout.flush()
        else:
            # get rid of any extra files now
            logs = os.listdir(os.getcwd());
            logs = [l for l in logs if '.log' in l]
            if len(logs)>0:
                for log in logs:
                    err = subprocess.call('rm %s'%log, shell=True)
             
    # since we moved dirs, going to put back to original directory
    os.chdir(original_dir)