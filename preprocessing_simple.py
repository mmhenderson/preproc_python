import numpy as np
import os, sys
import subprocess
import pandas as pd
from datetime import datetime
import copy
import scipy
import nibabel as nib

# change these depending on project, etc.
project_root = '/user_data/mmhender/data_UW/'
retino_path = '/lab_data/hendersonlab/retino_data/ANAT/'

# Code for preprocessing functional MRI data.
# Runs motion correction and de-trending

# SKIP registration step, this works if you have just one session.

# Based on matlab code from Rosanne Rademaker, adapted for python by MMH

# subject is your subject within individual experiment
# (a string like "S01")
# subject_FS is your subject in freesurfer/retinotopy directory.
# (a string like "AA")
# set debug=True if you're just testing the code, it will process first run then stop.

def run_preproc(subject, subject_FS, debug=False):
    
    debug=(int(debug)==1)
    print('subject=%s, subject_FS=%s, debug=%s'%(subject, subject_FS, debug))
    
    preproc_folder = os.path.join(project_root, 'DataPreprocSimple', subject)
    if not os.path.exists(preproc_folder):
        os.makedirs(preproc_folder)

    all_info_fn = os.path.join(preproc_folder, 'run_info_allsess.csv')
    if not os.path.exists(all_info_fn):
        get_session_info(subject)
        
    err = make_mc_template(subject, debug)
    if err: return
    
    # "MC"
    err = do_motion_correction(subject, debug)
    if err: return
    
    # "MC_DET"
    err = do_detrending(subject, debug)
    if err: return



def get_session_nums(subject):
    
    raw_folder = os.path.join(project_root, 'DataRaw', subject)
    sfolders = os.listdir(raw_folder)
    sess_nums = [int(s.split('Sess')[1]) for s in sfolders if 'Sess' in s]
    sess_nums = np.sort(np.array(sess_nums))
    
    return sess_nums



def get_session_info(subject):
    
    # this replaces the 'runs.list' in original matlab code. 
    
    preproc_folder = os.path.join(project_root, 'DataPreprocSimple', subject)
    if not os.path.exists(preproc_folder):
        os.makedirs(preproc_folder)

    all_info_fn = os.path.join(preproc_folder, 'run_info_allsess.csv')
    if os.path.exists(all_info_fn):
        print('Info file already exists at: %s'%(all_info_fn))
        print('If you want to make a new one, first delete this file.')
        return 0
    
    # Gather info about runs for each session.
    run_info_allsess = []

    sess_nums = get_session_nums(subject)

    for ss in sess_nums:

        raw_folder = os.path.join(project_root, 'DataRaw', subject, 'Sess%02d'%ss)

        # This is a pd dataframe that we already made.
        # It has info about which runs we did in which order.
        # Includes both field maps and functional scans.
        info_fn = os.path.join(raw_folder, 'Niftis', 'run_info.csv')
        print('loading info from %s'%info_fn)
        sys.stdout.flush()
        run_info = pd.read_csv(info_fn)

        # find just our functional runs here.
        func_inds = [i for i in np.arange(run_info.shape[0]) if \
                     'fmap' not in np.array(run_info['run_type'])[i]]

        # Creating a new DF that is across all sessions.
        rtmp = copy.deepcopy(run_info.iloc[func_inds])
        rtmp['num_in_session'] = np.arange(1, rtmp.shape[0]+1)
        rtmp['session'] = ss
        
        # figure out names of files that are input to subsequent processing
        raw_niis = [np.array(run_info['nifti_fn'])[rr] for rr in func_inds]

        # outstr = '_topup' # if we do unwarping, need this string
        outstr = '' # otherwise, empty string
        ext  = '.nii.gz'
        input_niis = np.array([r.split(ext)[0] + outstr + ext for r in raw_niis])

        rtmp['raw_fn'] = raw_niis
        rtmp['input_fn'] = input_niis
        
        run_info_allsess += [rtmp]

    run_info_allsess = pd.concat(run_info_allsess)
    run_info_allsess = run_info_allsess.set_index(np.arange(run_info_allsess.shape[0]))

    # allocating columns to store what my edited filenames will be.
    run_info_allsess['mc_fn'] = [' ' for rr in range(run_info_allsess.shape[0])]
    run_info_allsess['mc_det_fn'] = [' ' for rr in range(run_info_allsess.shape[0])]

    # Save this to a new file
    
    print(all_info_fn)
    run_info_allsess.to_csv(all_info_fn, index=False)
    
    
    
def make_mc_template(subject, debug=False):

    # We're going to take the first functional run from the first session, and grab
    # just the first timepoint of it. 
    # This is saved as our "motion correction template".

    preproc_folder = os.path.join(project_root, 'DataPreprocSimple', subject)
    
    mctemplate_file = os.path.join(preproc_folder, 'MCTemplateSess01.nii.gz')

   
    all_info_fn = os.path.join(preproc_folder, 'run_info_allsess.csv')
    run_info_allsess = pd.read_csv(all_info_fn)

    if os.path.exists(mctemplate_file):

        print('\nMotion correction template exists already, at: %s\n'%(mctemplate_file))

    else:

        print('\nMaking MC template\n')
        
        first_run_ind = np.where(np.array((run_info_allsess['session']==1) & \
                          (run_info_allsess['num_in_session']==1)))[0][0]

        first_nii_file = np.array(run_info_allsess['raw_fn'])[first_run_ind]

        # fslroi with (0, 1) means (start at time 0, take 1 volume)
        # https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Fslutils
        cmd = 'fslroi %s %s 0 1'%(first_nii_file, mctemplate_file)

        print(cmd)
        sys.stdout.flush()
        err = subprocess.call(cmd, shell=True)
        if err:
            print('Error in creating MC template.')
            return 0


def do_motion_correction(subject, debug=False):
    
    preproc_folder = os.path.join(project_root, 'DataPreprocSimple', subject)
    
    sess_nums = get_session_nums(subject)
    
    all_info_fn = os.path.join(preproc_folder, 'run_info_allsess.csv')
    run_info_allsess = pd.read_csv(all_info_fn)
    
    # this is our template that we're registering everything to.
    # important that it comes from the first session, so everything is ultimately in same space.
    mctemplate_file = os.path.join(preproc_folder, 'MCTemplateSess01.nii.gz')
        
    for ss in sess_nums:

        print('\nStarting motion correction for Sess%02d\n'%(ss))
        
        inds_do = np.where(np.array(run_info_allsess['session']==ss))[0]
        
        
        for ii in inds_do:

            rr = np.array(run_info_allsess['num_in_session'])[ii]
    
            if debug and ((ss>1) or (rr>1)):
                continue
                
            # going into my dataframe to get the filename
            # this was created in the last step, and now is input to this step.
            input_nii_file_full = np.array(run_info_allsess['input_fn'])[ii]
           
            # Naming the pre-processed files with a simpler naming scheme.
            # this makes it easier to find them later on.
            new_nii_file = 'Sess%02d_Run%02d_MC.nii.gz'%(ss, rr)
            new_nii_file_full = os.path.join(preproc_folder, new_nii_file)

            if os.path.exists(new_nii_file_full):

                print('\nMotion corrected volume for Run%02d exists already, at: %s'%(rr, new_nii_file_full))
                
                if run_info_allsess.at[ii,'mc_fn']!=new_nii_file:
                    # putting this new file name into my big dataframe, just to stay organized.
                    run_info_allsess.at[ii,'mc_fn'] = new_nii_file
                    # update the saved .csv file
                    print('writing to %s'%all_info_fn)
                    run_info_allsess.to_csv(all_info_fn, index=False)

            else:
    
                print('\nRunning motion correction for Run%02d...'%rr)
        
                # Use mcflirt to do motion correction.
                # using 6 dof here, but we could use 12 if desired.
                # https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=fsl;90247e02.1707
                cmd = 'mcflirt -in %s '%(input_nii_file_full) + \
                      '-smooth 0 -stages 4 -dof 6 ' + \
                      '-reffile %s '%(mctemplate_file) + \
                      '-out %s '%(new_nii_file_full) + \
                      '-plots -mats -v'

                print(cmd)
                sys.stdout.flush()
                err = subprocess.call(cmd, shell=True)

                if err:
                    print('Error with motion correction for %s'%(new_nii_file_full))
                    return 0
                else:
                    # putting this new file name into my big dataframe, just to stay organized.
                    run_info_allsess.at[ii,'mc_fn'] = new_nii_file
                    # update the saved .csv file
                    print('writing to %s'%all_info_fn)
                    run_info_allsess.to_csv(all_info_fn, index=False)

                    # Going to save out parameters for motion in this run.
                    # Later we can look at this for quality control.
                    # See: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/MCFLIRT
                    print('Gathering motion correction data (QC step)')
                    parfile = new_nii_file_full + '.par'
                    mat = read_motion_pars(parfile)

                    motion = pd.DataFrame(mat, columns=['x-rot (rad)', 'y-rot (rad)', 'z-rot (rad)', \
                                                        'x-trans (mm)', 'y-trans (mm)', 'z-trans (mm)'])

                    motiondir = os.path.join(preproc_folder, 'Motion')
                    if not os.path.exists(motiondir):
                        os.makedirs(motiondir)

                    motion_file = os.path.join(motiondir, 'Sess%02d_Run%02d_MC_motion.csv'%(ss, rr))
                    print(motion_file)
                    motion.to_csv(motion_file, index=False)

                    # going to gather first volume of each scan, after MC
                    # this will be used for quality control later on, by concatenating all of these.
                    firstvolsdir = os.path.join(preproc_folder, 'FirstVols')
                    if not os.path.exists(firstvolsdir):
                        os.makedirs(firstvolsdir)
                    firstvol_file = os.path.join(firstvolsdir, 'Sess%02d_Run%02d_MC_Vol1.nii.gz'%(ss, rr))

                    # take just first volume.
                    cmd = 'fslroi %s %s 0 1'%(new_nii_file_full, firstvol_file)
                    print(cmd)
                    err = subprocess.call(cmd, shell=True)

    firstvolsdir = os.path.join(preproc_folder, 'FirstVols')
    concat_file = os.path.join(firstvolsdir, 'concat_volume.nii.gz')
    files2concat = os.path.join(firstvolsdir, '*Vol1.nii.gz')
    print('Merging first volumes of all runs across all sessions (QC step)')
    cmd = 'fslmerge -t %s %s'%(concat_file, files2concat)
    print(cmd)
    err = subprocess.call(cmd, shell=True)

    
def read_motion_pars(par_file):

    with open(par_file) as f:
        lines = f.readlines()
        mat = []
        for line in lines:
            parts = line.split(' ')
            parts = [p for p in parts if (p!='') and (p!='\n')]
            nums = [float(p) for p in parts]
            mat += [nums]

    mat = np.array(mat)
    return mat

    
def do_detrending(subject, debug=False):
    
    preproc_folder = os.path.join(project_root, 'DataPreprocSimple', subject)
    
    sess_nums = get_session_nums(subject)
    
    all_info_fn = os.path.join(preproc_folder, 'run_info_allsess.csv')
    run_info_allsess = pd.read_csv(all_info_fn)
    
    for ss in sess_nums:

        print('\nStarting detrending for Sess%02d\n'%(ss))
        
        inds_do = np.where(np.array(run_info_allsess['session']==ss))[0]
        
        for ii in inds_do:

            rr = np.array(run_info_allsess['num_in_session'])[ii]
    
            if debug and ((ss>1) or (rr>1)):
                continue
                
            # going into my dataframe to get the filename
            # this was created in the last step, and now is input to this step.
            input_nii_file = np.array(run_info_allsess['mc_fn'])[ii]
            input_nii_file_full = os.path.join(preproc_folder, input_nii_file)
           
            # Naming the pre-processed files with a simpler naming scheme.
            # this makes it easier to find them later on.
            new_nii_file = 'Sess%02d_Run%02d_MC_DET.nii.gz'%(ss, rr)
            new_nii_file_full = os.path.join(preproc_folder, new_nii_file)

            if os.path.exists(new_nii_file_full):

                print('\nDetrended volume for Run%02d exists already, at: %s'%(rr, new_nii_file_full))
               
                if run_info_allsess.at[ii,'mc_det_fn']!=new_nii_file:
                    # putting this new file name into my big dataframe, just to stay organized.
                    run_info_allsess.at[ii,'mc_det_fn'] = new_nii_file
                    # update the saved .csv file
                    print('writing to %s'%all_info_fn)
                    run_info_allsess.to_csv(all_info_fn, index=False)

            else:
    
                print('\nRunning detrending for Run%02d...'%rr)
        
                high_pass_cutoff = 50; # volumes
                low_pass_cutoff = -1; # means actually no cutoff, it's just a highpass filter.

                # Note from RR code:
                # % Detrending (high pass filering) also removes the mean (the
                # % slowest drift in your data, if you will). If you want to add
                # % the mean back in you can use fslmaths to grab the mean first
                # % before you detrend (and add it back in during the detrend).
                mean_file_tmp = os.path.join(preproc_folder, 'tmpMean')
                cmd = 'fslmaths %s -Tmean %s'%(input_nii_file_full, mean_file_tmp)
                print(cmd)
                err = subprocess.call(cmd, shell=True)

                if err:
                    print('Error with detrending for %s'%(new_nii_file_full))
                    return 0
                
                # actual detrending
                # check out "fslmaths" in terminal, it prints info about all these options.
                # -bptf  <hp_sigma> <lp_sigma> : (-t in ip.c) Bandpass temporal filtering; 
                # nonlinear highpass and Gaussian linear lowpass (with sigmas in volumes, not seconds); 
                # set either sigma<0 to skip that filter
                cmd = 'fslmaths %s '%(input_nii_file_full) + \
                      '-bptf %d %d '%(high_pass_cutoff, low_pass_cutoff) + \
                      '-add %s '%(mean_file_tmp) + \
                      '%s '%(new_nii_file_full)

                print(cmd)
                err = subprocess.call(cmd, shell=True)

                if err:
                    print('Error with detrending for %s'%(new_nii_file_full))
                    return 0
                else:
                    # putting this new file name into my big dataframe, just to stay organized.
                    run_info_allsess.at[ii,'mc_det_fn'] = new_nii_file
                    # update the saved .csv file
                    print('writing to %s'%all_info_fn)
                    run_info_allsess.to_csv(all_info_fn, index=False)

