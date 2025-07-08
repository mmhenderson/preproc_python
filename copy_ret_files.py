import numpy as np
import os, sys
import subprocess

# change these depending on project, etc.
# project_root = '/user_data/mmhender/data_UW/'
# project_root = '/home/lab/hendersonlab/data_featsynth/'
project_root = '/lab_data/hendersonlab/data_featsynth/'
# retino_path = '/home/lab/hendersonlab/retino_data/'
retino_path = '/lab_data/hendersonlab/retino_data/'

# subject = 'S04'
# subject_FS = 'AE'
# debug=False
# sess=1

# copy_ret_files(subject, subject_FS, debug=False, sess=1)

def copy_all(subject, subject_FS, debug=False, sess=1):

    # helper functions for copying files from a project-specific folder 
    # into my main "retino_data" folder
    # change the paths above to set what those folders are
    # do this after unwarping, and before running any of the DOT_SCRIPTS code

    
    # subject is your subject within individual experiment
    # (a string like "S01")
    # subject_FS is your subject in freesurfer/retinotopy directory.
    # (a string like "AA")
    
    debug=(int(debug)==1)
    print('subject=%s, subject_FS=%s, debug=%s'%(subject, subject_FS, debug))
        
    err = copy_anat(subject, subject_FS, sess=sess)
    if err: return
    
    err = copy_functional(subject, subject_FS, sess=sess)
    if err: return
    
    err = copy_behav(subject, subject_FS, sess=sess)
    if err: return
    

def copy_anat(subject, subject_FS, sess=1):

    inpath = os.path.join(project_root, 'DataRaw', subject, 'Sess%d'%sess)
    dicom_dir = os.path.join(inpath, 'Dicoms')
    nifti_dir = os.path.join(inpath, 'Niftis')

    all_dirs = [d for d in os.listdir(dicom_dir)]
    all_dirs = [d for d in all_dirs if os.path.isdir(os.path.join(dicom_dir, d))]
    # all_dirs = [d.replace(' ', '_') for d in all_dirs]
    all_dirs = [d for d in all_dirs if ('anat-T1w_acq-MEMPRvNav_RMS' in d)]
    # ^ this is the averaged scan across echoes, which is what we want.

    new_dir = os.path.join(retino_path, 'RAW', subject_FS, 'NII')
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)
        
    for d in all_dirs:
        
        print(d)
        
        this_dir = os.path.join(dicom_dir, d)
        
        nii_files = os.listdir(this_dir)
        nii_files = [n for n in nii_files if '.nii.gz' in n and ' '  not in n]
        
        for n in nii_files:
            
            old_name_full = os.path.join(this_dir, n)
            new_name_full = os.path.join(new_dir, '%s_%s'%(d, n))

            # print('%s_%s'%(d, n))
            # print(new_name_full)
            # print(os.path.exists(new_name_full))
    
            if os.path.exists(new_name_full):
                print('%s already exists'%n)
            else:
                cmd = 'cp -p %s %s'%(old_name_full, new_name_full)
                print(cmd)
                err = subprocess.call(cmd, shell=True)

            

def copy_functional(subject, subject_FS, sess=1):

    inpath = os.path.join(project_root, 'DataRaw', subject, 'Sess%d'%sess)
    dicom_dir = os.path.join(inpath, 'Dicoms')
    nifti_dir = os.path.join(inpath, 'Niftis')

    all_dirs = [d for d in os.listdir(nifti_dir) if ' ' not in d]
    all_dirs = [d for d in all_dirs if os.path.isdir(os.path.join(nifti_dir, d))]
    # all_dirs = [d for d in all_dirs if ('vRF' in d and 'topup' not in d)]
    # copy over both raw and topuped data here, just to have it
    all_dirs = [d for d in all_dirs if (('vRF' in d))]

    new_dir = os.path.join(retino_path, 'RAW', subject_FS, 'NII')
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)

    
    for d in all_dirs:
        
        print(d)
        
        this_dir = os.path.join(nifti_dir, d)
        
        nii_files = os.listdir(this_dir)
        nii_files = [n for n in nii_files if '.nii.gz' in n and ' '  not in n and 'cal' not in n]
        
        for n in nii_files:
            
            old_name_full = os.path.join(this_dir, n)
            new_name_full = os.path.join(new_dir, '%s_%s'%(d, n))
    
            if os.path.exists(new_name_full):
                print('%s already exists'%n)
            else:
                cmd = 'cp -p %s %s'%(old_name_full, new_name_full)
                print(cmd)
                err = subprocess.call(cmd, shell=True)
    
    
                

def copy_behav(subject, subject_FS, sess=1):

    inpath = os.path.join(project_root, 'DataBehavior', subject, 'Sess%d'%sess)
    
    new_dir = os.path.join(retino_path, 'RAW', subject_FS, 'behav')
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)

    files = os.listdir(inpath)
    files = [f for f in files if 'RF_bar_multSize' in f]
    
    for f in files:

        old_name_full = os.path.join(inpath, f)
        new_name_full = os.path.join(new_dir, f)

        if os.path.exists(new_name_full):
                print('%s already exists'%new_name_full)
        else:
        
            cmd = 'cp -p %s %s'%(old_name_full, new_name_full)
            print(cmd)
            err = subprocess.call(cmd, shell=True)