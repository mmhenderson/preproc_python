import numpy as np
import os, sys
import pandas as pd
import nibabel as nib

# change this depending on project, etc.
project_root = '/home/lab/hendersonlab/data_featsynth/'
# project_root = '/lab_data/hendersonlab/data_featsynth/'


def get_data_concat_maintask(ss):

    # Concatenate the data from all runs and sessions of this task
    # Make a single file that has all runs, and save as .npy file
    # To make it smaller, we're going to mask out just the voxels in just the areas of interest.
    # (this is like MakeSamples in RR matlab code)
    
    ss = int(ss)
    
    subject = 'S%02d'%ss

    # hard code this. will need to adjust for different tasks.
    # this changes depending on experiment
    if ss<=2:     
        nTRs = 284 # first version
    else:
        nTRs = 344; # newer version

    
    # this is all the timing information. computed from behavioral .mat files.
    behav_data_folder = os.path.join(project_root, 'DataBehavior', subject)
    filename_load = os.path.join(behav_data_folder, '%s_maintask_behav_timing_info.csv'%subject)
    print('Reading from: %s'%filename_load)
    tr_info = pd.read_csv(filename_load)

    n_trs_total = tr_info.shape[0]
    
    # this is all the information about my scanning runs.
    preproc_folder = os.path.join(project_root, 'DataPreproc', subject)
    all_info_fn = os.path.join(preproc_folder, 'run_info_allsess.csv')
    run_info_allsess = pd.read_csv(all_info_fn)
    
    runs_maintask = run_info_allsess[run_info_allsess['run_type']=='vMain']

    # Now I'm sorting by run number within each session. 
    # Usually the runs are already in this order, but if we ever collected runs in the wrong order, 
    # like having to repeat a run at end of session, then we can fix it here. 
    # I always leave them in actual time order up to this point because it makes most sense with preprocessing 
    # procedure. But at this point we're done preprocessing, so put them in "run" order from now on...
    runs_maintask_sorted = runs_maintask.iloc[[]]
    sess = runs_maintask['date_raw']
    
    for ses in np.unique(sess):
    
        r = runs_maintask[sess==ses]
        sorder_runs = np.argsort(np.array(r['run_number']))
        if np.any(sorder_runs!=np.arange(len(sorder_runs))):
            print('warning: some runs were collected in different order. Switching them back now...')
            print(np.array(r['run_number']), sorder_runs)
        runs_maintask_sorted = pd.concat([runs_maintask_sorted, r.iloc[sorder_runs]])
    
    print(runs_maintask_sorted)
    

    # get ROI definitions here.
    roi_masks = get_ROI_masks(ss)
    roi_names = list(roi_masks.keys())
    hemis = list(roi_masks[roi_names[0]].keys())
    
    # the voxels we want to focus on are those that lie in any ROI of interest. 
    # will subdivide them more when we're actually doing analysis.
    roi_masks_concat = np.array([[roi_masks[r][h] for h in hemis] for r in roi_names])
    voxel_mask_all = np.any(np.any(roi_masks_concat, axis=0), axis=0)
    
    n_vox_total = len(voxel_mask_all)
    n_vox_use = np.sum(voxel_mask_all)
    
    preproc_data_folder = os.path.join(project_root, 'DataPreproc', subject)
    
    n_runs = len(runs_maintask)

    
    for ri in range(n_runs):
        
        # this is the final version of the file, after whatever steps of preprocessing have been done to it.
        final_preproc_file = os.path.join(preproc_data_folder, \
                                          np.array(runs_maintask_sorted['reg_mc_det_fn'])[ri])
        print('Loading from %s'%final_preproc_file)
        sys.stdout.flush()

        # try:
        dat = np.array(nib.load(final_preproc_file).get_fdata())
        assert(dat.shape[3]==nTRs)
        # except:
        #     continue
            
        # reshape the data. this has to match the ordering in the "ravel" command above.
        dat_reshape = np.reshape(dat, [n_vox_total, nTRs], order='C').T
        
        dat_use = dat_reshape[:, voxel_mask_all]
    
        print(dat_reshape.shape, dat_use.shape)
    
        if ri==0:
            dat_all = dat_use
        else:
            dat_all = np.concatenate([dat_all, dat_use], axis=0)
    
        print(dat_all.shape)

    assert(dat_all.shape[0]==n_trs_total)

    # save the array [nTRs x nVoxels] to a .npy file
    # calling this "DataConcat" because we have concatenated runs
    save_folder = os.path.join(project_root, 'DataConcat')
    if not os.path.exists(save_folder):
        os.makedirs(save_folder)
    save_data_fn = os.path.join(save_folder, '%s_maintask_concat.npy'%subject)

    print('Saving to %s'%save_data_fn)
    np.save(save_data_fn, {'dat_all': dat_all, \
                          'roi_masks': roi_masks, \
                          'voxel_mask_all': voxel_mask_all})


def get_ROI_masks(ss):

    subject = 'S%02d'%ss
        
    # gathering my ROI definitions and turning them into masks.
    
    roi_folder = os.path.join(project_root, 'DataPreproc', subject, 'ROI_labels')
        
    roi_names = np.array(['V1d','V1v','V2d','V2v','LOC'])
    hemis = np.array(['lh','rh'])
    
    roi_masks = dict([])
    
    for rname in roi_names:
        
        roi_masks[rname] = dict([])
        
        for hemi in hemis:
    
            mask_fn = os.path.join(roi_folder, '%s.%s.nii.gz'%(hemi, rname))
            
            mask = np.array(nib.load(mask_fn).get_fdata())
    
            # this is FLATTENED now - total number of voxels in one big list.
            roi_masks[rname][hemi] = (mask==1).ravel(order='C')
            
    
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
                    
                    print('%s-%s and %s-%s have %d voxels of overlap. Splitting them up evenly...'%\
                          (h, r1, h, r2, np.sum(overlap)))
    
                    orig_combined_size = np.sum(roi_masks[r1][h] | roi_masks[r2][h])
    
                    overlap_vox = np.where(overlap)[0]
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
                
                print('%s-%s and %s-%s have %d voxels of overlap. Splitting them up evenly...'%\
                      (h1, r1, h2, r2, np.sum(overlap)))

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