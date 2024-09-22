import numpy as np
import os, sys
import subprocess
import pandas as pd
from datetime import datetime
import copy
import scipy
import nibabel as nib

# change these depending on project, etc.
# project_root = '/user_data/mmhender/data_UW/'
project_root = '/lab_data/hendersonlab/data_featsynth/'
retino_path = '/lab_data/hendersonlab/retino_data/ANAT/'

# just make sure we're loading the right version of these modules
codepath = '/lab_data/hendersonlab/preproc_code/'
sys.path.insert(0, codepath)

from preproc_python import file_utils

# Code for preprocessing functional MRI data.
# Run cross-session alignment, motion correction, and de-trending

# Based on matlab code from Rosanne Rademaker, adapted for python and modified by MMH

# subject is your subject within individual experiment
# (a string like "S01")
# subject_FS is your subject in freesurfer/retinotopy directory.
# (a string like "AA")
# set debug=True if you're just testing the code, it will process first run then stop.

def run_preproc(subject, subject_FS, debug=False):
    
    debug=(int(debug)==1)
    print('subject=%s, subject_FS=%s, debug=%s'%(subject, subject_FS, debug))
    
    preproc_folder = os.path.join(project_root, 'DataPreproc', subject)
    if not os.path.exists(preproc_folder):
        os.makedirs(preproc_folder)

    all_info_fn = os.path.join(preproc_folder, 'run_info_allsess.csv')
    if not os.path.exists(all_info_fn):
        get_session_info(subject)

    # preparing registration files
    
    err = make_reg_templates(subject, debug)
    if err: return

    err = compute_sess2anat_reg(subject, subject_FS, debug)
    if err: return
        
    err = adjust_sess2anat_reg(subject, debug)
    if err: return

    err = make_mc_template_xfm(subject, debug)
    if err: return

    err = make_func2anat(subject, subject_FS)
    if err: return
    
    # below here are where we actually modify the data and start making preproc data files
    
    # niftis from this one has suffix "REG"
    err = do_registration(subject, debug)
    if err: return
    
    # "REG_MC"
    err = do_motion_correction(subject, debug)
    if err: return
    
    # "REG_MC_DET"
    err = do_detrending(subject, debug)
    if err: return



def get_session_nums(subject):
    
    raw_folder = os.path.join(project_root, 'DataRaw', subject)
    sfolders = os.listdir(raw_folder)
    sess_nums = [int(s.split('Sess')[1]) for s in sfolders if 'Sess' in s]
    sess_nums = np.sort(np.array(sess_nums))
    
    return sess_nums



def get_session_info(subject, do_topup = True):
    
    # this replaces the 'runs.list' in original matlab code.
    
    
    preproc_folder = os.path.join(project_root, 'DataPreproc', subject)
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

        # raw_folder = os.path.join(project_root, 'DataRaw', subject, 'Sess%02d'%ss)
        raw_folder = os.path.join(project_root, 'DataRaw', subject, 'Sess%01d'%ss)

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

        if do_topup:
            outstr = '_topup' # if we do unwarping, need this string
        else:
            outstr = '' # otherwise, empty string
        ext  = '.nii.gz'
        input_niis = np.array([r.split(ext)[0] + outstr + ext for r in raw_niis])

        if do_topup:
            # check the files exist (this is a check in case you forgot to run topups)
            for f in input_niis:
                if not os.path.exists(f):
                    raise RuntimeError('\n%s does not exist.\nMaybe you need to run unwarping first?\n'%f)
                
        rtmp['raw_fn'] = raw_niis
        rtmp['input_fn'] = input_niis
        
        run_info_allsess += [rtmp]

    run_info_allsess = pd.concat(run_info_allsess)
    run_info_allsess = run_info_allsess.set_index(np.arange(run_info_allsess.shape[0]))

    # allocating columns to store what my edited filenames will be.
    run_info_allsess['reg_fn'] = [' ' for rr in range(run_info_allsess.shape[0])]
    run_info_allsess['reg_mc_fn'] = [' ' for rr in range(run_info_allsess.shape[0])]
    run_info_allsess['reg_mc_det_fn'] = [' ' for rr in range(run_info_allsess.shape[0])]

    # Save this to a new file
    
    print(all_info_fn)
    run_info_allsess.to_csv(all_info_fn, index=False)
    
    
    
def make_reg_templates(subject, debug=False):

    preproc_folder = os.path.join(project_root, 'DataPreproc', subject)
    
    sess_nums = get_session_nums(subject)
    
    all_info_fn = os.path.join(preproc_folder, 'run_info_allsess.csv')
    run_info_allsess = pd.read_csv(all_info_fn)
    
    # Make templates for motion correction.
    for ss in sess_nums:
        
        if debug and ((ss>1)):
            continue

        # We're going to take the first functional run from this session, and grab
        # just the first timepoint of it. 
        # This is saved as our "registration template".
        # Used to align each session with our high-resolution anatomical scan.
        
        # Note these files are NOT exactly the template we'll use during mcflirt.
        # Later we make a file called MCTemplateXFMSess01, and that is the actual template to use.
        
        regtemplate_file = os.path.join(preproc_folder, 'RegTemplateSess%02d.nii.gz'%ss)

        if os.path.exists(regtemplate_file):
            
            print('\nRegistration template for Sess%02d exists already, at: %s\n'%(ss, regtemplate_file))

        else:
           
            print('\nMaking registration template for Sess%02d\n'%ss)
            first_run_ind = np.where(np.array((run_info_allsess['session']==ss) & \
                              (run_info_allsess['num_in_session']==1)))[0][0]

            first_nii_file = np.array(run_info_allsess['input_fn'])[first_run_ind]

            # fslroi with (0, 1) means (start at time 0, take 1 volume)
            # https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Fslutils
            cmd = 'fslroi %s %s 0 1'%(first_nii_file, regtemplate_file)
            
            print(cmd)
            sys.stdout.flush()
            err = subprocess.call(cmd, shell=True)
            if err:
                print('Error in creating registration template.')
                return 0
        
        
def compute_sess2anat_reg(subject, subject_FS, debug=False):

    # register each session to our high-res anatomy, using regtemplate
    
    preproc_folder = os.path.join(project_root, 'DataPreproc', subject)
    
    sess_nums = get_session_nums(subject)
    
    for ss in sess_nums:
        
        if debug and ((ss>1)):
            continue
        
        # I'm skipping the manual registration here - because we used auto-align during scanning, the volumes
        # are already quite close, and auto-registration should be sufficient. 
        # But make sure to double check it later!

        autoreg_file = os.path.join(preproc_folder, 'Sess%02d_raw2anat_xfm.dat'%ss)
        # autoreg_file = os.path.join(preproc_folder, 'registration_auto%02d.dat'%ss)
        # ^ this is the registration file we're going to save out. In "freesurfer" format
        # if you're loading the volumes into freeview, this registration would be used.

        bbout_file = os.path.join(preproc_folder, 'Sess%02d_raw2anat_bbout.nii.gz'%ss)
        # ^ this is the resampled template volume after registration.

        regtemplate_file = os.path.join(preproc_folder, 'RegTemplateSess%02d.nii.gz'%ss)
   
        if os.path.exists(autoreg_file):
        
            print('\nRegistration file for Sess%02d exists already, at: %s\n'%(ss, autoreg_file))

        else:
            
            print('\nRunning automatic registration (functional to anatomical) for Sess%02d\n'%ss)
            
            # calling freesurfer bbregister, boundary-based registration of volumes
            # https://surfer.nmr.mgh.harvard.edu/fswiki/bbregister
            # NOTE: skipping --init_reg flag here, bc we skipped manual step.
            
            # what we're specifically doing is aligning functional to the volume in:
            # ANAT/<subject_FS>/mri/orig.mgz
            # Freesurfer knows how to find this because you specify subject initials.
            
            cmd = 'bbregister --s %s '%(subject_FS) + \
                  '--mov %s '%(regtemplate_file) + \
                  '--reg %s '%(autoreg_file) + \
                  '--o %s '%(bbout_file) + \
                  '--bold'

            print(cmd)
            sys.stdout.flush()
            err = subprocess.call(cmd, shell=True)
            if err:
                print('Error in automatic boundary-based registration.')
                return 0

            # Take output from FreeSurfer's bbregister (.dat) and convert to FSL format (.mat)
            autoreg_file_fsl = os.path.join(preproc_folder, 'Sess%02d_raw2anat_xfm.mat'%ss)
            
            cmd = 'tkregister2 --s %s '%(subject_FS) + \
                  '--mov %s '%(regtemplate_file) + \
                  '--reg %s '%(autoreg_file) + \
                  '--fslregout %s '%(autoreg_file_fsl) + \
                  '--noedit'
            print(cmd)
            sys.stdout.flush()
            err = subprocess.call(cmd, shell=True)
            if err:
                print('Error in .dat to .mat conversion.')
                return 0

      
def adjust_sess2anat_reg(subject, debug=False):

    # This function will adjust the "raw2anat_xfm.mat" file so that it transforms
    # our data into a more reasonable space.
    # The raw version of this file is output by "bbregister", and it specifies a 
    # transform of our functional data (RegTemplate) to our high-res anatomy.
    # Because these volumes are different resolution and have different axes, 
    # the XFM swaps the axes. If we apply this XFM raw as the input to FLIRT, it 
    # makes it so the data goes out of our bounding box. So I'm editing the 
    # matrix so the output volumes are in similar coordinate space to originals.
    
    # In the original matlab code, this takes the place of RegHeaderCenterMOD workaround
    
    # Note that if scan protocols change, we might need to edit this.
    # Always check your volumes to make sure it's what you expect!
    # preproc_code/notebooks/check_alignment.ipynb
    
    preproc_folder = os.path.join(project_root, 'DataPreproc', subject)

    sess_nums = get_session_nums(subject)

    
    # This achieves a swapping of x/y columns, and negating of y column
    rot_mat = np.identity(4)[:,[0,2,1,3]]
    rot_mat[:,1] *= (-1)
    
    trans_rot_mat = copy.deepcopy(rot_mat)
    # this part may need to be edited someday
    # t = [-25, -190, -25]
    t = [-25, -25, 190] 
    # [R to L, I to S, A to P] ??
    trans_rot_mat[0:3,3] = t
    
    adj_mat_file = os.path.join(preproc_folder, 'xfm_adjustment.mat')
    
    print('Writing to %s'%adj_mat_file)
    file_utils.write_xfm_mat(trans_rot_mat, adj_mat_file)

    # this additional transformation MUST be the same for all sessions
    # otherwise they would no longer be aligned.


    for ss in sess_nums:

        if debug and ((ss>1)):
            continue
            
        autoreg_file = os.path.join(preproc_folder, 'Sess%02d_raw2anat_xfm.mat'%ss)
        xfm_file_new = os.path.join(preproc_folder, 'Sess%02d_raw2anat_xfm_ADJ.mat'%ss)

        if os.path.exists(xfm_file_new):
            
            print('\nFile exists already at: %s'%(xfm_file_new))
            print('Delete if you want to over-write it.')
            continue
            
        else:

            # concatenate the matrices with flirt
            # note this isn't just the same as multiplying them
            cmd = 'convert_xfm -omat %s '%(xfm_file_new) + \
                          '-concat %s %s '%(adj_mat_file, autoreg_file)
            print(cmd)
            sys.stdout.flush()
            err = subprocess.call(cmd, shell=True)

            if err:
                print('Error making adjusted matrix')
                return


def make_mc_template_xfm(subject, debug=False):

    preproc_folder = os.path.join(project_root, 'DataPreproc', subject)

    sess_nums = get_session_nums(subject)

    # reference, this affects header and how the file is loaded in FSLeyes and Freeview
    # make sure this is same for all sessions
    ref_file = os.path.join(preproc_folder, 'RegTemplateSess01.nii.gz')

    for ss in sess_nums:

        if debug and ((ss>1)):
            continue
            
        xfm_file_new = os.path.join(preproc_folder, 'Sess%02d_raw2anat_xfm_ADJ.mat'%ss)

        regtemplate_file = os.path.join(preproc_folder, 'RegTemplateSess%02d.nii.gz'%ss)

        # Make "MCTemplateXFM". 
        # This is the first volume of first run, transformed into the common space that
        # all our data eventually ends up in.
        # MCTemplateXFMSess01 is the target we use for motion correction with flirt.
        # If you load these files for all sessions, should be perfectly aligned.
        
        mctemplate_XFM_file = os.path.join(preproc_folder, 'MCTemplateXFMSess%02d.nii.gz'%ss)

        if os.path.exists(mctemplate_XFM_file):
            
            print('\nFile exists already at: %s'%(mctemplate_XFM_file))
            print('Delete if you want to over-write it.')
            continue
            
        else:
            
            cmd = 'flirt -in %s '%(regtemplate_file) + \
                '-ref %s '%(ref_file) + \
                '-applyxfm -init %s '%(xfm_file_new) + \
                '-out %s '%(mctemplate_XFM_file)
            print(cmd)
            sys.stdout.flush()
            err = subprocess.call(cmd, shell=True)
            if err:
                print('Error making MCTemplateXFM')
                return 0


def make_func2anat(subject, subject_FS):

    # making another registration
    # don't need this one until later when we're defining ROIs.
    # but making it now for convenience.
    
    preproc_folder = os.path.join(project_root, 'DataPreproc', subject)
        
    mctemplate_file = os.path.join(preproc_folder, 'MCTemplateXFMSess01.nii.gz')

    # create these files that describe the transformation
    # NOTE this is different from "raw2anat" file, because this one describes
    # transformation of the REG files (after first flirt transform is applied).
    autoreg_file = os.path.join(preproc_folder, 'Sess01_func2anat_xfm.dat')
    
    bbout_file = os.path.join(preproc_folder, 'Sess01_func2anat_bbout.nii.gz')
    # ^ this is the resampled template volume after registration.

    if os.path.exists(autoreg_file):
            
        print('\nFile exists already at: %s'%(autoreg_file))
        print('Delete if you want to over-write it.')

    else:
    
        cmd = 'bbregister --s %s '%(subject_FS) + \
                          '--mov %s '%(mctemplate_file) + \
                          '--reg %s '%(autoreg_file) + \
                          '--o %s '%(bbout_file) + \
                          '--bold'
        
        print(cmd)
        sys.stdout.flush()
        err = subprocess.call(cmd, shell=True)
        if err:
            print('Error creating func2anat file')
            return 0



def do_registration(subject, debug=False):
    
    preproc_folder = os.path.join(project_root, 'DataPreproc', subject)
    
    sess_nums = get_session_nums(subject)
    
    all_info_fn = os.path.join(preproc_folder, 'run_info_allsess.csv')
    run_info_allsess = pd.read_csv(all_info_fn)

    # reference, this affects header and how the file is loaded in FSLeyes and Freeview
    # make sure this is same for all sessions
    ref_file = os.path.join(preproc_folder, 'RegTemplateSess01.nii.gz')
        
    for ss in sess_nums:

        print('\nStarting coregistration (applying transformation) for Sess%02d\n'%(ss))
        
        inds_do = np.where(np.array(run_info_allsess['session']==ss))[0]
        
        # same transformation matrix gets applied for all runs in this session.
        xfm_file = os.path.join(preproc_folder, 'Sess%02d_raw2anat_xfm_ADJ.mat'%ss)

        for ii in inds_do:

            rr = np.array(run_info_allsess['num_in_session'])[ii]
            
            if debug and ((ss>1) or (rr>1)):
                continue
                
            # going into my dataframe to get the filename
            # this was created in the last step, and now it's the input to this step.
            input_nii_file = np.array(run_info_allsess['input_fn'])[ii]
            
            # Naming the pre-processed files with a simpler naming scheme.
            # this makes it easier to find them later on.
            new_nii_file = 'Sess%02d_Run%02d_REG.nii.gz'%(ss, rr)
            new_nii_file_full = os.path.join(preproc_folder, new_nii_file)

            if os.path.exists(new_nii_file_full):

                print('\nCoregistered volume for Run%02d exists already, at: %s'%(rr, new_nii_file_full))
                
                if run_info_allsess.at[ii,'reg_fn']!=new_nii_file:
                    # putting this new file name into my big dataframe, just to stay organized.
                    run_info_allsess.at[ii,'reg_fn'] = new_nii_file
                    # update the saved .csv file
                    print('writing to %s'%all_info_fn)
                    run_info_allsess.to_csv(all_info_fn, index=False)

            else:
    
                print('\nRunning coregistration for Run%02d...'%rr)
        
                # Using flirt to do the transformation to this subject's native space.
                # this just ensures all sessions are in same space, before we start motion correction.
                cmd = 'flirt -in %s '%(input_nii_file) + \
                      '-ref %s '%(ref_file) + \
                      '-out %s '%(new_nii_file_full) + \
                      '-init %s '%(xfm_file) + \
                      '-applyxfm -v'
                print(cmd)
                sys.stdout.flush()
                err = subprocess.call(cmd, shell=True)

                if err:
                    print('Error with transformation for %s'%(new_nii_file_full))
                    return 0
                else:
                    # putting this new file name into my big dataframe, just to stay organized.
                    run_info_allsess.at[ii,'reg_fn'] = new_nii_file
                    # update the saved .csv file
                    print('writing to %s'%all_info_fn)
                    run_info_allsess.to_csv(all_info_fn, index=False)

                    
def do_motion_correction(subject, debug=False):
    
    preproc_folder = os.path.join(project_root, 'DataPreproc', subject)
    
    sess_nums = get_session_nums(subject)
    
    all_info_fn = os.path.join(preproc_folder, 'run_info_allsess.csv')
    run_info_allsess = pd.read_csv(all_info_fn)
    
    # this is our template that we're registering everything to.
    # important that it comes from the first session, so everything is ultimately in same space.
    mctemplate_file = os.path.join(preproc_folder, 'MCTemplateXFMSess01.nii.gz')
        
    for ss in sess_nums:

        print('\nStarting motion correction for Sess%02d\n'%(ss))
        
        inds_do = np.where(np.array(run_info_allsess['session']==ss))[0]
        
        
        for ii in inds_do:

            rr = np.array(run_info_allsess['num_in_session'])[ii]
    
            if debug and ((ss>1) or (rr>1)):
                continue
                
            # going into my dataframe to get the filename
            # this was created in the last step, and now is input to this step.
            input_nii_file = np.array(run_info_allsess['reg_fn'])[ii]
            input_nii_file_full = os.path.join(preproc_folder, input_nii_file)
           
            # Add another suffix here
            new_nii_file = 'Sess%02d_Run%02d_REG_MC.nii.gz'%(ss, rr)
            new_nii_file_full = os.path.join(preproc_folder, new_nii_file)

            if os.path.exists(new_nii_file_full):

                print('\nMotion corrected volume for Run%02d exists already, at: %s'%(rr, new_nii_file_full))
                
                if run_info_allsess.at[ii,'reg_mc_fn']!=new_nii_file:
                    # putting this new file name into my big dataframe, just to stay organized.
                    run_info_allsess.at[ii,'reg_mc_fn'] = new_nii_file
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
                    run_info_allsess.at[ii,'reg_mc_fn'] = new_nii_file
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

                    motion_file = os.path.join(motiondir, 'Sess%02d_Run%02d_REG_MC_motion.csv'%(ss, rr))
                    print(motion_file)
                    motion.to_csv(motion_file, index=False)

                    # going to gather first volume of each scan, after MC
                    # this will be used for quality control later on, by concatenating all of these.
                    firstvolsdir = os.path.join(preproc_folder, 'FirstVols')
                    if not os.path.exists(firstvolsdir):
                        os.makedirs(firstvolsdir)
                    firstvol_file = os.path.join(firstvolsdir, 'Sess%02d_Run%02d_REG_MC_Vol1.nii.gz'%(ss, rr))

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
    
    preproc_folder = os.path.join(project_root, 'DataPreproc', subject)
    
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
            input_nii_file = np.array(run_info_allsess['reg_mc_fn'])[ii]
            input_nii_file_full = os.path.join(preproc_folder, input_nii_file)
           
            # Add another suffix here
            new_nii_file = 'Sess%02d_Run%02d_REG_MC_DET.nii.gz'%(ss, rr)
            new_nii_file_full = os.path.join(preproc_folder, new_nii_file)

            if os.path.exists(new_nii_file_full):

                print('\nDetrended volume for Run%02d exists already, at: %s'%(rr, new_nii_file_full))
               
                if run_info_allsess.at[ii,'reg_mc_det_fn']!=new_nii_file:
                    # putting this new file name into my big dataframe, just to stay organized.
                    run_info_allsess.at[ii,'reg_mc_det_fn'] = new_nii_file
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
                    run_info_allsess.at[ii,'reg_mc_det_fn'] = new_nii_file
                    # update the saved .csv file
                    print('writing to %s'%all_info_fn)
                    run_info_allsess.to_csv(all_info_fn, index=False)

