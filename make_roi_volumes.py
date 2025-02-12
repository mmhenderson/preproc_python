import nibabel as nib
import numpy as np
import os, sys
import subprocess
import scipy.io

project_root = '/home/lab/hendersonlab/data_featsynth/'
retino_path = '/home/lab/hendersonlab/retino_data/'

def proc_retino_annot(subject_FS):

    print(subject_FS)

    for hemi in ['lh', 'rh']:
            
        # Annot file is starting point - file that has all the retino areas defined in it.
        # This file should live in:
        #<retino_path>/ANAT/<subject_FS>/label
        annot_filename = '%s.%s_retinotopy'%(hemi, subject_FS)
        
        # New files go here.
        # I'm making a new folder "labels" to keep things together
        outdir = os.path.join(retino_path, 'ANAT', subject_FS, 'labels')
        if not os.path.exists(outdir):
            os.makedirs(outdir)
            
        cmd = 'mri_annotation2label '+\
              '--subject %s '%(subject_FS)+\
              '--hemi %s '%(hemi)+\
              '--annotation %s '%(annot_filename) +\
              '--outdir %s '%(outdir)
                
        print(cmd)
        sys.stdout.flush()
        err = subprocess.call(cmd, shell=True)

def proc_floc_annot(subject_FS):
    
    print(subject_FS)

    for hemi in ['lh', 'rh']:
            
        # Annot file is starting point - file that has all the fLoc areas defined in it.
        # This file should live in:
        #<retino_path>/ANAT/<subject_FS>/label
        annot_filename = '%s.%s_floc'%(hemi, subject_FS)
        
        # New files go here.
        # I'm making a new folder "labels" to keep things together
        outdir = os.path.join(retino_path, 'ANAT', subject_FS, 'labels')
        if not os.path.exists(outdir):
            os.makedirs(outdir)
            
        cmd = 'mri_annotation2label '+\
              '--subject %s '%(subject_FS)+\
              '--hemi %s '%(hemi)+\
              '--annotation %s '%(annot_filename) +\
              '--outdir %s '%(outdir)
                
        print(cmd)
        sys.stdout.flush()
        err = subprocess.call(cmd, shell=True)

    
def make_vols_func(subject, subject_FS):

    print(subject, subject_FS)

    roi_names = ['V1d','V1v','V2d','V2v', 'LOC']
        
    # this folder has defs of my labels in anatomical/surface space
    # drawn in freeview
    labels_folder_anat = os.path.join(retino_path, 'ANAT', subject_FS, 'labels')

    # this folder is where i'm putting my new nifti files
    # which are in volume space, for this particular project
    labels_folder = os.path.join(project_root, 'DataPreproc', subject, 'ROI_labels')
    if not os.path.exists(labels_folder):
        os.makedirs(labels_folder)
    
    
    for rname in roi_names:
            
        for hemi in ['lh', 'rh']:
                       
            labelfile = os.path.join(labels_folder_anat, '%s.%s.label'%(hemi, rname))
            
            # template file: this is the motion correction template for very first session.
            # this defines the space my functional data from this project are in.
            templatefile = os.path.join(project_root, 'DataPreproc', subject, 'MCTemplateXFMSess01.nii.gz')
            
            
            # registration file: this is the file that defines a mapping from my functional space
            # (for this subject, in this project) to my anatomical space (high-res anatomical scan
            # for this subject). we made this file during preprocessing.py
            regfile = os.path.join(project_root, 'DataPreproc', subject, 'Sess01_func2anat_xfm.dat')
            
            out_volume = os.path.join(labels_folder, '%s.%s.nii.gz'%(hemi, rname))
            
            fillthresh = 0.30;
            
            cmd = 'mri_label2vol '+\
                  '--subject %s '%(subject_FS) +\
                  '--hemi %s '%(hemi) +\
                  '--label %s '%(labelfile) +\
                  '--reg %s '%(regfile) +\
                  '--temp %s '%(templatefile) +\
                  '--o %s '%(out_volume) +\
                  '--fillthresh %.1f '%(fillthresh) +\
                  '--proj frac 0 1 0.1 '
                    
            print(cmd)
            sys.stdout.flush()
            err = subprocess.call(cmd, shell=True)

            