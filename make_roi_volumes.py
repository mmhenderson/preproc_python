import os, sys

print('\npython version = %s'%(sys.version))
print('env is: %s'%sys.prefix)

import nibabel as nib
import numpy as np

import subprocess
import scipy.io

from scipy.stats import t

project_root = '/lab_data/hendersonlab/data_featsynth/'
retino_path = '/lab_data/hendersonlab/retino_data/'
# project_root = '/home/lab/hendersonlab/data_featsynth/'
# retino_path = '/home/lab/hendersonlab/retino_data/'

codepath = '/lab_data/hendersonlab/code_featsynth/preproc_code/'
sys.path.insert(0, codepath)

from preproc_python import analyze_categ_loc


def proc_retino_annot(subject_FS):

    # as of april 2025, not using this. i am making the .label files one at a time.

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

    # as of april 2025, not using this. i am making the .label files one at a time.
    
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

    roi_names = ['V1v','V1d','V2v','V2d','V3v','V3d', 'V3A','V3B','hV4',\
                 'IPS0','IPS1','IPS2','IPS3','IPS4','LO1','LO2','VO1','VO2', \
                'LO', 'pFus']
        
    # this folder has defs of my labels in anatomical/surface space
    # drawn in freeview
    # labels_folder_anat = os.path.join(retino_path, 'ANAT', subject_FS, 'labels')
    labels_folder_retino = os.path.join(retino_path, 'ANAT', subject_FS, 'label', 'retino')
    labels_folder_floc = os.path.join(retino_path, 'ANAT', subject_FS, 'label', 'floc')
    # if you save the labels to a different folder, need to change this.

    # this folder is where i'm putting my new nifti files
    # which are in volume space, for this particular project
    labels_folder_new = os.path.join(project_root, 'DataPreproc', subject, 'ROI_labels')
    if not os.path.exists(labels_folder_new):
        os.makedirs(labels_folder_new)
    
    
    for rname in roi_names:
            
        for hemi in ['lh', 'rh']:

            if rname in ['LO', 'pFus']:
                labelfile = os.path.join(labels_folder_floc, '%s.%s.label'%(hemi, rname))
            else:
                labelfile = os.path.join(labels_folder_retino, '%s.%s.label'%(hemi, rname))
            
            # template file: this is the motion correction template for very first session.
            # this defines the space my functional data from this project are in.
            templatefile = os.path.join(project_root, 'DataPreproc', subject, 'MCTemplateXFMSess01.nii.gz')
            
            
            # registration file: this is the file that defines a mapping from my functional space
            # (for this subject, in this project) to my anatomical space (high-res anatomical scan
            # for this subject). we made this file during preprocessing.py
            regfile = os.path.join(project_root, 'DataPreproc', subject, 'Sess01_func2anat_xfm.dat')
            
            out_volume = os.path.join(labels_folder_new, '%s.%s.nii.gz'%(hemi, rname))
            
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


def threshold_ROI_masks(subject):
    
    # gathering my ROI definitions and turning them into masks.
    
    roi_folder = os.path.join(project_root, 'DataPreproc', subject, 'ROI_labels')

    roi_names = ['V1v','V1d','V2v','V2d','V3v','V3d', 'V3A','V3B','hV4',\
                 'IPS0','IPS1','IPS2','IPS3','IPS4','LO1','LO2','VO1','VO2', \
                'LO', 'pFus']
    
    # This is a list of which contrast you will use to define which areas.
    # FFA is defined based on face > all, and so on.
    # You must double check this if it's a new experiment!! 
    # this comes from the FEAT analysis, there is a .png file that has the design.

    contrast_dict = analyze_categ_loc.get_contrast_names(subject)
   
    hemis = np.array(['lh','rh'])
    
    for rname in roi_names:

        if rname in ['LO','pFus']:
            # thresholding the object-selective LOC area. going to combine these subregions
            cope_name = contrast_dict['object_gr_scrambled']
        else:
            # The fLoc stimuli are the same size as my standard main task stimuli. 10 deg.
            # So we can use this task additionally as a way of thresholding the retinotopic ROIs, for responsive voxels.
            cope_name = contrast_dict['stim_gr_baseline']

        # print(rname, cope_name)
            
        for hemi in hemis:

            # this is the original ROI. Directly from what we drew on the brain.
            # it may include non-significant voxels.
            mask_fn = os.path.join(roi_folder, '%s.%s.nii.gz'%(hemi, rname))
            # print(mask_fn)
            roi_orig = nib.load(mask_fn)
            
            mask = np.array(roi_orig.get_fdata())
    
            print('\nBefore thresh, %s-%s has %d voxels'%(hemi, rname, np.sum(mask)))

            
            # this is where we do the significance test, get t-threshold.
            vox_above_thresh, tstats, my_t = analyze_categ_loc.get_tstats_masked(subject, cope_name)
            
            roi_thresh = (mask==1) & vox_above_thresh
    
            print('After thresh, %s-%s has %d voxels'%(hemi, rname, np.sum(roi_thresh)))
            

            # make my new ROI file - just like the old one but fewer voxels.
            roi_new = nib.Nifti1Image(roi_thresh, roi_orig.affine, roi_orig.header)
            
            mask_fn_new = os.path.join(roi_folder, '%s.%s_THRESH.nii.gz'%(hemi, rname))
            print(mask_fn_new)
            nib.save(roi_new, mask_fn_new)




def get_ROI_masks(ss, threshold = True, fix_overlap = True, remove_overlapping = True):

    subject = 'S%02d'%ss
        
    # gathering my ROI definitions and turning them into masks.
    
    roi_folder = os.path.join(project_root, 'DataPreproc', subject, 'ROI_labels')

    # choose how to combine the sub-parts here, and any you want to ignore.
    roi_subparts = {'V1': ['V1v','V1d'], \
                'V2': ['V2v','V2d'], \
                'V3': ['V3v','V3d'], \
                'hV4': ['hV4'], \
                'V3AB': ['V3A','V3B'], \
                'IPS': ['IPS0','IPS1','IPS2','IPS3','IPS4'], \
                'LOC': ['LO','pFus']}
    
    roi_names = list(roi_subparts.keys())
    
    hemis = np.array(['lh','rh'])
    
    roi_masks = dict([])
    
    for rname in roi_names:
        
        roi_masks[rname] = dict([])

        rname_list  = roi_subparts[rname]
        
        for hemi in hemis:

            mask_list = []
            
            for r in rname_list:
                    
                if threshold:
                    mask_fn = os.path.join(roi_folder, '%s.%s_THRESH.nii.gz'%(hemi, r))
                else:
                    mask_fn = os.path.join(roi_folder, '%s.%s.nii.gz'%(hemi, r))

                # print(mask_fn)
                mask = np.array(nib.load(mask_fn).get_fdata())
                # print(np.sum(mask))
                
                mask_list += [mask]

            mask = np.stack(mask_list, axis=0)
            
            mask = np.any(mask==1, axis=0)

            # this is FLATTENED now - total number of voxels in one big list.
            roi_masks[rname][hemi] = (mask).ravel(order='C')
            # print('%s-%s: total vox %d'%(rname, hemi, np.sum(mask)))
            
    if fix_overlap:
            
        # now i'm going to handle the overlap between ROIs.
        # this happens because ROIs are defined in surface space, and they get projected back into 
        # volume space. this can lead to overlap between neighboring areas, even if they didn't look overlapping 
        # in our surface space.
        # generally will just split up voxels evenly (half and half, randomly) when there is overlap.
        # if we're using category selective areas (like LOC, FFA, and such), and also retinotopic areas, they 
        # can overlap a lot. So we might have to make a different choice about how to 
        # deal with overlap between them. 
        
        # first dealing within the overlap within-hemifield.
        
        for h in hemis:
            
            for ri1, r1 in enumerate(roi_names):
                
                for r2 in roi_names[ri1+1:]:
        
                    overlap = roi_masks[r1][h] & roi_masks[r2][h]
                    
                    if np.sum(overlap)>0:
                        
                        print('%s: %s-%s and %s-%s have %d voxels of overlap.'%\
                              (subject, h, r1, h, r2, np.sum(overlap)))
        
                        orig_combined_size = np.sum(roi_masks[r1][h] | roi_masks[r2][h])
        
                        overlap_vox = np.where(overlap)[0]

                        if remove_overlapping:

                            print('removing any overlapping vox...')
                            
                            # this is the strictest thing you can do, remove anything that is overlapping.
                            # for fLoc, this might be best (??)
                            roi_masks[r1][h][overlap] = False
                            roi_masks[r2][h][overlap] = False

                            # make sure this worked and we really fixed the overlap
                            overlap = roi_masks[r1][h] & roi_masks[r2][h]
                            assert(np.sum(overlap)==0)

                            # make sure nothing else weird happened
                            new_combined_size = np.sum(roi_masks[r1][h] | roi_masks[r2][h])
                            print(orig_combined_size, new_combined_size, (orig_combined_size - len(overlap_vox)))
                            
                            assert(new_combined_size==(orig_combined_size - len(overlap_vox)))
                            
            
                        elif (r1=='LOC') or (r2=='LOC'):
                            
                            print('putting all the voxels in the retinotopic ROI...')
                                  
                            if r1=='LOC':
                                overlap1 = []
                                overlap2 = overlap_vox
                            else:
                                overlap1 = overlap_vox
                                overlap2 = []
                            print(len(overlap1), len(overlap2))
                            
                            # remove them from the opposite roi
                            roi_masks[r1][h][overlap2] = False
                            roi_masks[r2][h][overlap1] = False
            
                            # make sure this worked and we really fixed the overlap
                            overlap = roi_masks[r1][h] & roi_masks[r2][h]
                            assert(np.sum(overlap)==0)
            
                            # make sure nothing else weird happened
                            new_combined_size = np.sum(roi_masks[r1][h] | roi_masks[r2][h])
                            assert(new_combined_size==orig_combined_size)
                                
                        else:

                            print('splitting them up evenly...')
                            
                            # less strict - keep them but only in one ROI.
                            # splitting every other vox. this is good because it's reproducible.
                            # you could also do it randomly, but you would have to save the random seed. 
                            # very important that we know exactly which voxels we grabbed here.
                            overlap1 = overlap_vox[0::2]
                            overlap2 = overlap_vox[1::2]
                            print(len(overlap1), len(overlap2))
                
                            # remove them from the opposite roi
                            roi_masks[r1][h][overlap2] = False
                            roi_masks[r2][h][overlap1] = False
            
                            # make sure this worked and we really fixed the overlap
                            overlap = roi_masks[r1][h] & roi_masks[r2][h]
                            assert(np.sum(overlap)==0)
            
                            # make sure nothing else weird happened
                            new_combined_size = np.sum(roi_masks[r1][h] | roi_masks[r2][h])
                            assert(new_combined_size==orig_combined_size)
            
    
    # now look for overlap across-hemifields.
    # testing all possible pairs of ROIs for overlap.
    # note that the majority of these have very little chance at overlapping.
    # it's really just the bordering ones. 
    
    h1 = hemis[0]
    h2 = hemis[1]
    
    for r1 in roi_names:
        
        for r2 in roi_names:
    
            overlap = roi_masks[r1][h1] & roi_masks[r2][h2]
    
            if np.sum(overlap)>0:
                
                print('%s: %s-%s and %s-%s have %d voxels of overlap. Splitting them up evenly...'%\
                      (subject, h1, r1, h2, r2, np.sum(overlap)))

                orig_combined_size = np.sum(roi_masks[r1][h1] | roi_masks[r2][h2])
    
                overlap_vox = np.where(overlap)[0]
                # splitting every other voxel again
                overlap1 = overlap_vox[0::2]
                overlap2 = overlap_vox[1::2]
                print(len(overlap1), len(overlap2))
    
                # remove them from the opposite roi
                roi_masks[r1][h1][overlap2] = False
                roi_masks[r2][h2][overlap1] = False
                
                overlap = roi_masks[r1][h1] & roi_masks[r2][h2]
                assert(np.sum(overlap)==0)

                # make sure nothing else weird happened
                new_combined_size = np.sum(roi_masks[r1][h1] | roi_masks[r2][h2])
                assert(new_combined_size==orig_combined_size)

    
    return roi_masks