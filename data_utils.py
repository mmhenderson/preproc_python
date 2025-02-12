import numpy as np
import os, sys
import pandas as pd
import nibabel as nib
import scipy.stats

# change this depending on project, etc.
project_root = '/home/lab/hendersonlab/data_featsynth/'
# project_root = '/lab_data/hendersonlab/data_featsynth/'


def load_main_task_labels(ss):
    
    behav_fn = os.path.join(project_root, 'DataBehavior', 'S%02d'%ss, \
                                          'S%02d_maintask_behav_trial_info.csv'%ss)
    print('loading from %s'%behav_fn)
    bdat = pd.read_csv(behav_fn)

    return bdat


def load_main_task_data(ss, make_time_resolved = True):
    
    """
    Load trial-by-trial data for all main task trials for a single subject (ss)
    make_time_resolved is bool for whether to return the tr-by-tr signal in 
    addition to averaged signal for each trial. 
    
    returns main_data (dictionary with ROIs as keys, holds all the data)
    
    """

    subject = 'S%02d'%ss
        
    load_folder = os.path.join(project_root, 'DataConcat')
    load_data_fn = os.path.join(load_folder, '%s_maintask_concat.npy'%subject)
    print('Loading data from: %s'%(load_data_fn))
    d = np.load(load_data_fn, allow_pickle=True).item()
    n_vox_all = d['dat_all'].shape[1]
    raw_data = d['dat_all']


    # this is all the timing information. computed from behavioral .mat files.
    behav_data_folder = os.path.join(project_root, 'DataBehavior', subject)
    filename_load = os.path.join(behav_data_folder, '%s_maintask_behav_timing_info.csv'%subject)
    print('Reading from: %s'%filename_load)
    tr_info = pd.read_csv(filename_load)
    
    n_trs_total = tr_info.shape[0]
    
    # info about individual trials
    filename_load = os.path.join(behav_data_folder, '%s_maintask_behav_trial_info.csv'%subject)
    print('Reading from: %s'%filename_load)
    trial_info = pd.read_csv(filename_load)
    
    n_trials_total = trial_info.shape[0]

    # defining some timing stuff...
    
    # how many TRs in the task
    if ss<=2:     
        nTRs = 284 # first version
    else:
        nTRs = 344; # newer version

    
    
    # which TRs am i averaging over? From the target onset time.
    avgTRs_targ = [3, 7];
    # in this expt: TRs are 1.0 seconds in length.
    # TRs go like...TR 0 = 0.0 seconds (time at which stim onset happened), 
    # TR 1 = 1.0 seconds after stim onset, and so on
    
    # if we're returning time-resolved data, how many TRs to include per trial?
    # nTRs_concat = 7;
    nTRs_concat = 14;

    # First: we will zscore the data from each run to normalize it. 
    # this corrects for any big differences in signals across runs
    # Note that after detrending, we put the mean back in. so there can be big differences 
    # in mean as well as variance.
    zdata = np.full(fill_value = np.nan, shape=raw_data.shape)
    
    run_labels = tr_info['run_num_overall']
    print(np.unique(run_labels))
    for run in np.unique(run_labels):
    
        run_inds = run_labels==run
        assert(np.sum(run_inds)==nTRs)
    
        zdata[run_inds,:] = scipy.stats.zscore(raw_data[run_inds,:], axis=0)

    assert(not np.any(np.isnan(zdata)))
    
    # Next i'm going to "epoch" the data by finding which TRs go with which trials
    # gather data for each individual trial.
    # average over a fixed time window following stim onset.
    dat_avg = np.full(fill_value=np.nan, shape=(n_trials_total, n_vox_all))
    
    if make_time_resolved:
        dat_by_tr = np.full(fill_value=np.nan, shape=(n_trials_total, nTRs_concat, n_vox_all))
    else:
        dat_by_tr = None

    # find onset of each trial
    # it's where the "stim_on" variable goes from 0 to 1. 
    diffs = np.diff(np.array([0.0] + list(tr_info['stim_on'])))
    stim_onsets = diffs==1
    assert(np.all(np.array(tr_info['stim_on'])[stim_onsets]==1))
    assert(np.sum(stim_onsets)==n_trials_total)

    tc = -1 # count over trials 
    
    # we're going to do this one run at a time, even though the data is concatenated.
    # this is important because it makes sure we never accidentally combine data from two runs.
    # like for example if you're looking at an event close to the end of run 1, you might accidentally 
    # grab TRs from start of run 2 in the concatenated data. but that would be wrong.
    # so this helps keep the runs separate.
    for ri, rr in enumerate(np.unique(tr_info['run_num_overall'])):
    
        run_inds = tr_info['run_num_overall']==rr
        assert(np.sum(run_inds)==nTRs)
    
        cur_dat = zdata[run_inds]
    
        these_stim_onsets = np.where(stim_onsets[run_inds])[0]
    
        # print(np.max(these_stim_onsets))
        for tt in these_stim_onsets:
    
            tc+=1
            
            # this is the time window following the stimulus onset
            # idea is this should roughly capture peak of bold signal. 
            # you can verify this by plotting the data timecourses
            window_start = tt + avgTRs_targ[0]
            window_stop = tt + avgTRs_targ[1]
    
            assert(window_stop<nTRs)
    
            # take the average of timepoints within my time window.
            # average across time, within each voxel separately
            dat_avg[tc, :] = np.mean(cur_dat[window_start:window_stop, :], axis=0)
    
            if make_time_resolved:

                for tr in range(nTRs_concat):
                            
                    # if this gives an error, it means your nTRs_concat window might be too long
                    # the window has to be short enough that it is defined even for the last trial
                    # assert((tt+nTRs_concat) < nTRs)
                    if tt + tr >= nTRs:
                        print('Warning: TR %d is past end of the run, for trial %d'%(tr, tt))
                        # we're putting nans in here.
                        # not good!! don't do this really
                        dat_by_tr[tc, tr, :] = np.nan
                        
                    else:

                        # taking the whole epoch of data corresponding to this trial
                        dat_by_tr[tc, tr, :] = cur_dat[tt + tr, :]
        
    assert(tc==(n_trials_total-1))

    # Now that we have the "epoched" data, for all voxels, we're going to mask out the ROIs of interest.
    # Make a big dictionary where keys are ROIs. 
    main_data = dict([])
    
    # these are the areas we want to look at
    roi_names = ['V1','V2','LOC']
    hemis = ['lh', 'rh']
    
    for rname in roi_names:
    
        # the actual key names in the roi_masks dict also include some
        # subsets of ROIs, like d=dorsal and v=ventral.
        # we are going to combine those sub-parts here.
        # also combine lh/rh here.
    
        keys_include = [r for r in list(d['roi_masks'].keys()) if rname in r]
    
        print('%s, including:'%rname)
        print(keys_include)
        masks = []
    
        for kk in keys_include:
            for hh in hemis:
    
                # the roi mask is defined over all voxels, in whole brain
                mask_big = d['roi_masks'][kk][hh]
                # we want just the part of this that corresponds to the data
                mask_small = mask_big[d['voxel_mask_all']]
    
                print('    %s-%s has %d voxels'%(kk, hh, np.sum(mask_small)))
                masks += [mask_small]
    
        masks_concat = np.array(masks).astype(int)
    
        # checking for overlap between sub-parts of the ROI
        assert(not np.any(np.sum(masks_concat, axis=0)>1))

        # want voxels that are defined in any sub-part
        mask_use = np.any(masks_concat, axis=0)
    
        print('%s has %d total voxels'%(rname, np.sum(mask_use)))
    
        # finally, use the ROI mask to pick my voxels from the epoched data.
        main_data[rname] = dict([])
        main_data[rname]['dat_avg'] = dat_avg[:, mask_use]
        main_data[rname]['dat_by_tr'] = dat_by_tr[:, :, mask_use]

    
    return main_data
