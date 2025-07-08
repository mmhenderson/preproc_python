import numpy as np
import os, sys
import pandas as pd
from datetime import datetime
import copy

# change this depending on project, etc.
# project_root = '/home/lab/hendersonlab/data_featsynth/'
project_root = '/lab_data/hendersonlab/data_featsynth/'

# make sure we're loading the right version of these modules
# codepath = '/home/lab/hendersonlab/code_featsynth/preproc_code/'
codepath = '/lab_data/hendersonlab/code_featsynth/preproc_code/'
sys.path.insert(0, codepath)

from preproc_python import file_utils


def proc(ss):

    ss = int(ss)
    # ss is your subject number (like 1, or 2, etc)
    
    # this makes:
    # 'maintask_behav_run_info.csv'
    process_runs_maintask(ss)

    # this checks the above file against the scanner data files
    cross_check_behav_fmri_maintask(ss)

    # this makes:
    # 'maintask_behav_trial_info.csv'
    process_trials_maintask(ss)

    # this makes:
    # 'maintask_behav_timing_info.csv'
    process_timing_maintask(ss)
    

def process_runs_maintask(ss):

    # this is where we gather info about the individual runs in the experiment.
    
    # find all main task behavioral (.mat) files, and make a dataframe with info about the files.
    # not actually processing details of each run yet, just filenames/when they were collected.
    
    # save a csv file n_runs long
    
    # subject should be like: 1, 2 etc

    subject = 'S%02d'%ss
    
    # where is the behav data?
    behav_data_folder = os.path.join(project_root, 'DataBehavior', subject)
    
    # these are the session folders I have.
    # Note that these session folder numbers do not necessarily align with the session in task. 
    # BC retinotopy is usually Sess1 folder, so task starts in Sess2.
    sfolders = os.listdir(behav_data_folder)
    sess_nums = [int(s.split('Sess')[1]) for s in sfolders if 'Sess' in s]
    sess_nums = np.sort(np.array(sess_nums))
    
    # name of the task that we're analyzing.
    # this is a string that has to be in all the filenames of your behavioral .mat files
    # can change this to match your filenames for your task.
    task_name = 'maintask_fMRI_pilot'
    
    behav_run_info = pd.DataFrame(columns=['file_name', 'sess_in_task', 'run_in_task',\
                                           'num_in_session', 'sess_actual', \
                                         'date_raw','dstring',\
                                         'time_raw', 'tstring', \
                                         'length_sec', 'avg_acc'])
    
    n_sessions_expected = 3;
    n_runs_expected = 12;
    
    all_files = []
    all_num_in_sess = []
    all_sess_actual = []
    
    # this is a loop over whichever physical sessions were performed (incl retinotopy)
    for ses in sess_nums:
        
        subfolder = os.path.join(behav_data_folder, 'Sess%d'%(ses))
        files = os.listdir(subfolder) 
        # print(files)
        # print(task_name)
        files = [f for f in files if task_name in f]
        # print(files)
        print('Sess %d: found %d behav files for main task'%(ses, len(files)))
    
        if (ses>1) & (len(files)!=n_runs_expected):
            print('WARNING: num runs is different than expected for %s session %d. Check this!'%(subject, ses))
    
        all_files += [files]
        
        time_raw = [f.split('.mat')[0][-6:] for f in files]
        time_raw = np.array([float(t) for t in time_raw])
    
        # apply argsort twice, it tells us the rank of each element
        torder = np.argsort(np.argsort(time_raw)) + 1
        print(time_raw, torder)
    
        all_num_in_sess += [torder]
        all_sess_actual += [ses * np.ones(len(files),)]
        
    
    all_files = np.concatenate(all_files)
    all_num_in_sess = np.concatenate(all_num_in_sess)
    all_sess_actual = np.concatenate(all_sess_actual)

    ii = 0;
    # now this is a loop over "session in task", which always goes 1,2,3
    for sess_in_task in np.arange(1, n_sessions_expected+1):
    
        for run_in_task in np.arange(1, n_runs_expected+1):
    
            file_index = np.where(['Sess%02d_Run%02d'%(sess_in_task, run_in_task) in file for file in all_files])[0]
            if len(file_index)==0:
                print('WARNING: Sess%02d Run%02d is missing, skipping this'%(sess_in_task, run_in_task))
                continue
            
            # There has to be exactly one of these runs. If there are two, this means you have a duplicate.
            # Usually this would be if there was a mistake with one, so you need to remove that one.
            assert(len(file_index)==1)
                
            file_index = file_index[0]
            
            f = all_files[file_index]
    
            print(f)
            
            date_raw = f.split('.mat')[0][-13:-7]
            dstring = datetime.strftime(datetime.strptime(date_raw, '%y%m%d'), '%D')
            
            time_raw = f.split('.mat')[0][-6:]
            tstring = datetime.strftime(datetime.strptime(time_raw, '%H%M%S'), '%I:%M %p')
            
            subject_check = f.split('_')[0]
            assert(subject_check==subject)
    
            subfolder = os.path.join(behav_data_folder, 'Sess%d'%(all_sess_actual[file_index]))
            
            print(os.path.join(subfolder, f))
            
            p = file_utils.load_mat_behav_data(os.path.join(subfolder, f), varname='p')
            assert(len(p)==1)
            p = p[0]
            
            ii+=1
            # sess_in_task and run_in_task are both related to what was entered in the experimental task code.
            # num_in_session and sess_actual are related to the actual date on which runs were collected.
            # So, for some subjects Sess01 is actually the retinotopy session, and Sess02 is the first 
            # session of the main task. Sess02 would have: sess_actual=2, sess_in_task=1
            # And sometimes we might have had to split functional runs for a given "session" 
            # across days - for example collecting runs 1-13 in session 3, and runs 14-16 in session 4. 
            # then sess_actual would be different for those days, but sess_in_task would be the same.
            behav_run_info = pd.concat([behav_run_info, pd.DataFrame({'file_name': f, \
                                                     'sess_in_task': sess_in_task, \
                                                     'run_in_task': run_in_task,\
                                                     'num_in_session': all_num_in_sess[file_index], \
                                                     'sess_actual': all_sess_actual[file_index], \
                                                     'date_raw': date_raw, \
                                                     'dstring': dstring, \
                                                     'time_raw': time_raw, \
                                                     'tstring': tstring, \
                                                     'length_sec': p['total_exp_time'], \
                                                     'avg_acc': p['accuracy'], \
                                                        }, \
                                                       index = [ii])])
        
    print('All sessions: %d runs of main task total'%(behav_run_info.shape[0]))

    filename_save = os.path.join(behav_data_folder, '%s_maintask_behav_run_info.csv'%subject)
    print('Saving to: %s'%filename_save)
    behav_run_info.to_csv(filename_save, index=False)


# def process_runs_maintask(ss):

#     # this is where we gather info about the individual runs in the experiment.
    
#     # find all main task behavioral (.mat) files, and make a dataframe with info about the files.
#     # not actually processing details of each run yet, just filenames/when they were collected.
    
#     # save a csv file n_runs long
    
#     # subject should be like: 1, 2 etc

#     subject = 'S%02d'%ss
    
#     # where is the behav data?
#     behav_data_folder = os.path.join(project_root, 'DataBehavior', subject)
    
#     sfolders = os.listdir(behav_data_folder)
#     sess_nums = [int(s.split('Sess')[1]) for s in sfolders if 'Sess' in s]
#     sess_nums = np.sort(np.array(sess_nums))
    
#     # name of the task that we're analyzing.
#     # this is a string that has to be in all the filenames of your behavioral .mat files
#     # can change this to match your filenames for your task.
#     task_name = 'maintask_fMRI_pilot'
    
#     behav_run_info = pd.DataFrame(columns=['file_name', 'sess_in_task', 'run_in_task',\
#                                            'num_in_session', 'sess_actual', \
#                                          'date_raw','dstring',\
#                                          'time_raw', 'tstring', \
#                                          'length_sec', 'avg_acc'])

#     n_runs_expected = 12;
    
#     for ses in sess_nums:
        
#         subfolder = os.path.join(behav_data_folder, 'Sess%d'%(ses))
#         files = os.listdir(subfolder) 
#         # print(files)
#         # print(task_name)
#         files = [f for f in files if task_name in f]
#         # print(files)
#         print('Sess %d: found %d behav files for main task'%(ses, len(files)))

#         if len(files)==0:
#             continue
            
#         date_raw = [f.split('.mat')[0][-13:-7] for f in files]
#         dstring = [datetime.strftime(datetime.strptime(d, '%y%m%d'), '%D') for d in date_raw]
    
        
#         time_raw = [f.split('.mat')[0][-6:] for f in files]
#         tstring = [datetime.strftime(datetime.strptime(t, '%H%M%S'), '%I:%M %p') for t in time_raw]
        
#         sess_in_task = [int(f.split('Sess')[1][0:2]) for f in files]
    
#         run_in_task = [int(f.split('Run')[1][0:2]) for f in files]
#         if not np.all(np.sort(run_in_task)==np.arange(1, n_runs_expected+1)):
#             print('WARNING: you may have a missing or duplicated file for %s session %s'%(subject, ses))
#             print(np.sort(run_in_task))
#         # assert(np.all(np.sort(run_in_task)==np.arange(1, n_runs_expected+1)))
#         # make sure there is exactly one of each run number
#         # if this condition is not met, we need to check the behav data to see what's wrong.
#         # it might be that we missed a run, or duplicated a run somehow.

        
#         subject_check = np.array([f.split('_')[0] for f in files])
#         # print(subject_check, subject)
#         assert(np.all(subject_check==subject))
    
#         # sort by time at which they were collected
#         sorder_time = np.argsort(np.array(time_raw))
       
#         # the run numbers within a given session should always be increasing by steps of 1.
#         # check that this is true. if not, could mean we collected runs in wrong order. 
#         # This actually happens sometimes if we needed to re-do a run, so it can be fine.
#         # We just have to make sure that the mri data and behav data are aligned right, which
#         # will be addressed later in preproc code.
#         if not np.all(np.diff(np.array(run_in_task)[sorder_time])==1):
#             print('WARNING: you may have either a duplicate or transposed file for %s session %s'%(subject, ses))
#             print(np.array(run_in_task)[sorder_time])

#         # now sorting by ascending run number
#         sorder = np.argsort(run_in_task)

        
#         for ii, fi in enumerate(sorder):

#             print(os.path.join(subfolder, files[fi]))
            
#             p = file_utils.load_mat_behav_data(os.path.join(subfolder, files[fi]), varname='p')
#             assert(len(p)==1)
#             p = p[0]
    
#             # sess_in_task and run_in_task are both related to what was entered in the experimental task code.
#             # num_in_session and sess_actual are related to the actual date on which runs were collected.
#             # So, for some subjects Sess01 is actually the retinotopy session, and Sess02 is the first 
#             # session of the main task. Sess02 would have: sess_actual=2, sess_in_task=1
#             # And sometimes we might have had to split functional runs for a given "session" 
#             # across days - for example collecting runs 1-13 in session 3, and runs 14-16 in session 4. 
#             # then sess_actual would be different for those days, but sess_in_task would be the same.
#             behav_run_info = pd.concat([behav_run_info, pd.DataFrame({'file_name': files[fi], \
#                                                      'sess_in_task': sess_in_task[fi], \
#                                                      'run_in_task': run_in_task[fi],\
#                                                      'num_in_session': ii+1, \
#                                                      'sess_actual': ses, \
#                                                      'date_raw': date_raw[fi], \
#                                                      'dstring': dstring[fi], \
#                                                      'time_raw': time_raw[fi], \
#                                                      'tstring': tstring[fi], \
#                                                      'length_sec': p['total_exp_time'], \
#                                                      'avg_acc': p['accuracy'], \
#                                                         }, \
#                                                        index = [ii])])
    
    
#     print('All sessions: %d runs of main task total'%(behav_run_info.shape[0]))

#     filename_save = os.path.join(behav_data_folder, '%s_maintask_behav_run_info.csv'%subject)
#     print('Saving to: %s'%filename_save)
#     behav_run_info.to_csv(filename_save, index=False)

def cross_check_behav_fmri_maintask(ss):

    subject = 'S%02d'%ss
    # check the "maintask_behav_run_info.csv" that I just made, against my "run_info_allsess.csv" file
    # make sure the runs line up right.
    
    # Load the info about behavior files. this is a CSV file,  generated based on what is in my DataBehavior folders.
    behav_data_folder = os.path.join(project_root, 'DataBehavior', subject)
    filename_load = os.path.join(behav_data_folder, '%s_maintask_behav_run_info.csv'%subject)
    behav_run_info = pd.read_csv(filename_load)

    # sort the array now by TIME of when files were collected
    # this is usually same as its usual order, EXCEPT if the files were ever collected in a funny order (like if 
    # we have to redo a run at end of session). The above .csv file always goes in run number order ascending, 
    # while the MRI data file .csv file (loaded next) goes in actual time order. 
    # so we want to put this in actual time order so we can compare against mri data files.
    behav_run_info_timesorted = behav_run_info.iloc[[]]
    date = behav_run_info['date_raw']
    
    for d in np.unique(date):
    
        print(d)
        b = behav_run_info[date==d]
        sorder_time = np.argsort(b['time_raw'])
        behav_run_info_timesorted = pd.concat([behav_run_info_timesorted, b.iloc[sorder_time]])

        
    # now I'm going to load info about my MRI files.
    # these are in a separate CSV file that we made as part of preprocessing pipeline.
    # want to make sure we have correspondence between each main task file in behavior and MRI dataframes.
    preproc_folder = os.path.join(project_root, 'DataPreproc', subject)
    all_info_fn = os.path.join(preproc_folder, 'run_info_allsess.csv')
    run_info_allsess = pd.read_csv(all_info_fn)

    # main task only
    m = run_info_allsess[run_info_allsess['run_type'] == 'vMain']
    
    # check number of files
    assert(m.shape[0]==behav_run_info_timesorted.shape[0])
    
    # check dates of files
    np.all(np.array(m['dstring'])==np.array(behav_run_info_timesorted['dstring']))
    np.all(np.array(m['session'])==np.array(behav_run_info_timesorted['sess_actual']))

    # checking that the lengths line up
    tr_length = 1.0;
    sec_check = np.array(m['nTRs']) * tr_length
    assert(np.all(np.abs(sec_check - np.array(behav_run_info_timesorted['length_sec']))<0.50))
    
    # check that the run numbers line up.
    # for MRI data: run number is the actual number in the filename, which is entered during the scanning session.
    # for behavior data: this is the number that we enter into the matlab script.
    # if these are not perfectly lined up, might indicate we have a duplicate, or entered the wrong run number somewhere.
    print(np.array(m['run_number']))
    print(np.array(behav_run_info_timesorted['run_in_task']))
    
    assert(np.all(np.array(m['run_number'])==np.array(behav_run_info_timesorted['run_in_task'])))
    
    # check times... this is a bit more tricky
    # How many seconds were the start times offset by?
    # this is really just a feature of the exact setup at BRIDGE - i think the clocks are different
    # between linux computer where behavior is run from, and the MRI computer.
    # it ends up looking like the behavior file starts slightly later (not actually true because you
    # start behav script first). but clocks are offset. so it starts about 2 min after.
    # this is just a sanity check that we're always using same behavior file that corresponds to same mri file.
    # as long as the same offset range holds, this will work. but for another setup, likely will have to edit this.
    t_mri = [datetime.strptime('%s'%t, '%H%M%S') for t in np.array(m['time_raw'])]
    t_behav = [datetime.strptime('%s'%t, '%H%M%S') for t in np.array(behav_run_info_timesorted['time_raw'])]
    
    sec_offset = np.array([(tb - tm).seconds for [tm, tb] in zip(t_mri, t_behav)])
    print('behav and MRI time offsets in sec are:')
    # NOTE sometimes when there are longer offsets here, this can be because we gave subject a longer "break" time
    # in this case, scanner may have created file before matlab was launched?
    print(sec_offset)
    # this cutoff is kind of arbitrary, if it errors out just check it
    assert(np.all(sec_offset>0) and np.all(sec_offset<60*6))
    # assert(np.all(sec_offset>0) and np.all(sec_offset<60*4))

    print('Done, check successful!')


def process_trials_maintask(ss):

    # this is where we will gather information about every trial in the experiment
    # save a .csv file that is n_trials long
    
    # ss = 2;
    subject = 'S%02d'%ss
    
    acc_thresh = 0.52
    
    # Load the info about behavior files. this is a CSV file,  generated based on what is in my DataBehavior folders.
    behav_data_folder = os.path.join(project_root, 'DataBehavior', subject)
    filename_load = os.path.join(behav_data_folder, '%s_maintask_behav_run_info.csv'%subject)
    behav_run_info = pd.read_csv(filename_load)
    
    # going to copy all the keys directly over from my run_info dataframe
    # this is stuff that is the same for a whole run
    run_info_keys = list(behav_run_info.keys())
    
    # these are keys that i will directly copy out of my saved matlab behavioral data
    # they correspond to elements in the "p" structure, which are lists n_trials long
    keys_get = ['miniblock_number', 'miniblock_task', 'correct_resp', 'image_names', \
                'super_names', 'basic_names', 'left_names','right_names', \
                'cue_level', 'imfns', 'response', 'resp_time_from_onset']
    
    # full list of dataframe columns to create
    allkeys = ['subject', 'ss', 'trial_num_in_run', 'trial_num_overall', 'run_num_overall', \
               'correct', 'ex_num', 'image_type_name', \
               'correct_name', 'distract_name'] + \
            keys_get + run_info_keys
    bdat = pd.DataFrame(columns = allkeys)
    
    tc = 0
    # loop over all the runs we collected
    # for ii in range(2):
    for ii in range(behav_run_info.shape[0]):
        
        ses = np.array(behav_run_info['sess_actual'])[ii]
        subfolder = os.path.join(behav_data_folder, 'Sess%d'%(ses))
        
        fn = os.path.join(subfolder, np.array(behav_run_info['file_name'])[ii])
        
        print('Loading from %s'%(fn))
        
        p = file_utils.load_mat_behav_data(fn, varname='p')
        assert(len(p)==1)
        p = p[0]

        if ss<=2:
            n_runs = 16;
        else:
            n_runs = 12;
        assert((np.mod(p['run_num_overall']-1, n_runs)+1)==np.array(behav_run_info['run_in_task'])[ii])
        
        acc_check = np.mean(np.array(p['response'])==np.array(p['correct_resp']))
        
        assert(p['accuracy']==acc_check)
    
                
        no_resp_trials = (p['resp_time_from_onset']==0) | np.isnan(p['resp_time_from_onset'])
    
        rr = p['run_num_overall']
        print('Run %d: Accuracy is: %.2f'%(rr, acc_check))
        print('Run %d: Number of no-response trials is: %d'%(rr, np.sum(no_resp_trials)))
    
        if acc_check<acc_thresh:
            # we could choose to ignore bad runs like this. for now, keep them but print a warning
            print('Warning: accuracy for run %d is below threshold of %.2f!'%(rr, acc_thresh))
            
        # checking some more stuff about names, responses
        names_lr = np.array([np.array(p['left_names']), np.array(p['right_names'])]).T
        
        super_task_inds = np.array(p['cue_level'])=='super'
        basic_task_inds = np.array(p['cue_level'])=='basic'
        
        # what was the exact name of "correct" response category?
        sname = np.array(p['super_names'])
        bname = np.array(p['basic_names'])
        correct_name = [b if bi else s for [s,b,bi] in zip(sname, bname, basic_task_inds)]
        
        # correct_name = np.array(p['super_names'])
        # correct_name[basic_task_inds] = np.array(p['basic_names'])[basic_task_inds]
        
        # check that "correct_resp" field makes sense
        # 1=left, 2=right
        # print(p['basic_names'])
        # print(np.array(p['basic_names']))
        # print(np.array(p['basic_names'])[basic_task_inds])
        # print(correct_name[basic_task_inds])
        
        # print(correct_name, names_lr)
        lr_correct = [np.where(c==n)[0][0] for c, n in zip(correct_name, names_lr)]
        assert(np.all(np.array(lr_correct)+1 == np.array(p['correct_resp'])))
        
        # what was name of "incorrect" (i.e. distractor) category?
        distract_name = np.array([n[(c!=n)][0] for c, n in zip(correct_name, names_lr)])
        assert(not np.any(correct_name==distract_name))
        
        p['num_trials']
        n_trials = p['num_trials']
    
        # now loop over individual trials in the run
        for ti in range(n_trials):
    
            tc+=1
            
            # going to make a new dataframe for this trial/this row.
            # start with the elements that are already defined in the run_info df
            newdf = copy.deepcopy(behav_run_info.iloc[ii:ii+1])
            
            # now adding a bunch more elements
            newdf.insert(1, 'subject', subject)
            newdf.insert(1, 'ss', ss)
            newdf.insert(1, 'run_num_overall', p['run_num_overall'])
            newdf.insert(1, 'trial_num_overall', tc)
            newdf.insert(1, 'trial_num_in_run', ti+1)
            
            for kk in keys_get:
            
                assert(len(np.array(p[kk]))==n_trials)
                newdf.insert(1, kk, np.array(p[kk])[ti])
            
            correct = np.array(p['response'])[ti]==np.array(p['correct_resp'])[ti]
            newdf.insert(1, 'correct', correct)
            
            # more specific stuff about the image shown here
            bname = np.array(p['image_names'])[ti].split('/')[0]
            assert(bname==np.array(p['basic_names'])[ti])
            
            ex_num = int(np.array(p['image_names'])[ti].split('/')[1][2:])
            newdf.insert(1, 'ex_num', ex_num)
            
            image_type_name = np.array(p['image_names'])[ti].split('/')[2].split('.png')[0]
            newdf.insert(1, 'image_type_name', image_type_name)
        
            newdf.insert(1, 'correct_name', correct_name[ti])
            newdf.insert(1, 'distract_name', distract_name[ti])
        
            bdat = pd.concat([bdat, newdf])
    
    filename_save = os.path.join(behav_data_folder, '%s_maintask_behav_trial_info.csv'%subject)
    print('Saving to: %s'%filename_save)
    bdat.to_csv(filename_save, index=False)


def process_timing_maintask(ss):

    # this is where we process information about individual TR timing in the experiment
    # save a csv file that is nTRs long
    # (this is longer than the number of trials)
    
    subject = 'S%02d'%ss

    # this changes depending on experiment
    if ss<=2:     
        nTRs = 284 # first version
    else:
        nTRs = 344; # newer version
    tr_dur = 1.0 # this too
    
    # Find onset times for all my TRs
    TR_onset = np.arange(0, nTRs*tr_dur, tr_dur)
    TR_middles = TR_onset + tr_dur/2 # this marks the "center" of each TR in seconds
        
    # Load the info about behavior files
    behav_data_folder = os.path.join(project_root, 'DataBehavior', subject)
    filename_load = os.path.join(behav_data_folder, '%s_maintask_behav_run_info.csv'%subject)
    behav_run_info = pd.read_csv(filename_load)
    
    trial_count = 0
    
    # ii = 0;
    # for ii in range(2):
    for ii in range(behav_run_info.shape[0]):
        
        ses = np.array(behav_run_info['sess_actual'])[ii]
        subfolder = os.path.join(behav_data_folder, 'Sess%d'%(ses))
        
        fn = os.path.join(subfolder, np.array(behav_run_info['file_name'])[ii])
        
        print('Loading from %s'%(fn))
        
        p = file_utils.load_mat_behav_data(fn, varname='p')
        assert(len(p)==1)
        p = p[0]
        
        # gathering info about timing of events within this run.
        flip_times = np.array(p['stim_flips']) - p['start_exp_time']
        n_trials = flip_times.shape[0]
        # print(n_trials)
        
        time_tol = 0.10
        
        # checking that the times make sense. at beginning, we have some events that are fixed length
        # first trial should start after this.
        expected_first_trial = p['start_fix_sec'] + p['miniblock_instr_sec'] + p['miniblock_delay_sec']
        assert(np.abs(flip_times[0,0] - expected_first_trial)<time_tol)
        
        # make sure events lasted right length
        event_lengths = np.diff(flip_times, axis=1)
        
        assert(np.all(np.abs(event_lengths[:,0]-p['stim_time_sec'])<time_tol))
        assert(np.all(np.abs(event_lengths[:,1]-p['delay_time_sec'])<time_tol))
        assert(np.all(np.abs(event_lengths[:,2]-p['cue_time_sec'])<time_tol))
        
        expected_end_time = flip_times[n_trials-1,:][3] + p['iti_sec'][n_trials-1] + p['end_fix_sec']
        assert(np.abs(p['total_exp_time'] - expected_end_time)<time_tol)
        
        # now I'm going to make lists that are nTRs long
        # describe what happens on each individual trial in experiment.
        # we'll need this to align the stimulus sequence and MRI data.
        
        trial_num_overall = np.zeros([nTRs,], dtype=int)
        stim_on = np.zeros([nTRs,], dtype=bool)
        cue_on = np.zeros([nTRs,], dtype=bool)
        
        
        for ti in range(n_trials):
            
            stim_onset = flip_times[ti,0]
            # which TR is this a part of? look for whichever TR has its "middle" closest to my event.
            stim_tr_ind = np.argmin(np.abs(stim_onset - TR_middles))
        
            stim_offset = flip_times[ti,1]
            stimoff_tr_ind = np.argmin(np.abs(stim_offset - TR_middles))
        
            if stim_tr_ind==stimoff_tr_ind:
                # this can happen because the events are shorter than 1 TR. 
                # so we'll mark the event as lasting 1 TR long. otherwise we would miss it.
                stimoff_tr_ind = stim_tr_ind + 1
        
            # print(stim_tr_ind, stimoff_tr_ind)
            # now mark which TRs are part of this event
            stim_on[stim_tr_ind:stimoff_tr_ind] = True
            
            cue_onset = flip_times[ti,2]
            cue_tr_ind = np.argmin(np.abs(cue_onset - TR_middles))
        
            cue_offset = flip_times[ti,3]
            cueoff_tr_ind = np.argmin(np.abs(cue_offset - TR_middles))
        
            if cue_tr_ind==cueoff_tr_ind:
                # this shouldn't actually happen
                cueoff_tr_ind = cue_tr_ind + 1
        
            # now mark which TRs are part of this event
            cue_on[cue_tr_ind:cueoff_tr_ind] = True
        
            # now i'm labeling all of these TRs as belonging to the current trial number.
            # it lasts from stim onset to cue offset (doesn't count the ITI)
            # will need this later on
            trial_count +=1
            trial_num_overall[stim_tr_ind:cueoff_tr_ind] = trial_count
        
        
        # this will be true for events that last 1 TR or less.
        # if events are longer, there would be more TRs marked.
        # print(np.sum(stim_on), n_trials, trial_count)
        assert(np.sum(stim_on)==n_trials)
    
        # putting it together
        run_num_overall = np.ones([nTRs,]) * (ii+1) # this just lists the run number in case we need it
        tr_info_this_run = pd.DataFrame(np.array([stim_on, cue_on, trial_num_overall, run_num_overall, \
                                                  TR_onset, TR_middles]).T, \
                                        columns = ['stim_on','cue_on','trial_num_overall','run_num_overall', \
                                                   'TR_onset', 'TR_middles'])
        
        if ii==0:
            tr_info = tr_info_this_run
        else:
            tr_info = pd.concat([tr_info, tr_info_this_run])
    
        
    filename_save = os.path.join(behav_data_folder, '%s_maintask_behav_timing_info.csv'%subject)
    print('Saving to: %s'%filename_save)
    tr_info.to_csv(filename_save, index=False)


def get_categ_info():

    # Creating a file that contains info about the 48 categories used in the experiment.
    # It's based on a file that originally had info about 64 categories, and we reduced it to 48 
    # for the fMRI study. Just adjusting the file so we can load it easily.
    # Should only need to run this once, same for all subjects.
    # (except for S02, which was the pilot subject w 64 categories).

    ecoset_info_path = '/user_data/mmhender/stimuli/ecoset_info/'

    fn = os.path.join(ecoset_info_path, 'categ_use_ecoset.npy')
    info = np.load(fn, allow_pickle=True).item()
    bnames = np.array(list(info['binfo'].keys()))
    snames = np.array(list(info['sinfo'].keys()))

    basic_use_6 = {'insect': ['beetle','bee','butterfly','grasshopper','caterpillar','moth'],
             'mammal': ['dog', 'squirrel', 'elephant', 'cow', 'pig', 'rabbit'],
             'vegetable': ['pea', 'corn', 'onion', 'cabbage', 'beet', 'asparagus'],
             'fruit': ['grape', 'cherry', 'raspberry', 'pear', 'banana', 'coconut'],
             'tool': ['pencil', 'knife', 'broom', 'hammer', 'shovel', 'scissors'],
             'musical instrument': ['bell', 'piano', 'violin', 'trumpet', 'clarinet', 'cymbal'],
             'furniture': ['table', 'bench', 'couch', 'television', 'bed', 'lamp'],
             'vehicle': ['train', 'airplane', 'car', 'bus', 'motorcycle', 'canoe']}

    assert( np.all( snames==np.array(list(basic_use_6.keys())) ) )

    info_48 = dict([])
    info_48['sinfo'] = dict([])
    info_48['binfo'] = dict([])

    for si in range(len(snames)):

        info_48['sinfo'][snames[si]] = {'super_name': snames[si], \
                                       'basic_names': basic_use_6[snames[si]]}

        for bname in basic_use_6[snames[si]]:

            info_48['binfo'][bname] = info['binfo'][bname]

    fn2save = os.path.join(codepath, 'preproc_python', 'categ_info_48.npy')
    print(fn2save)
    np.save(fn2save, info_48)