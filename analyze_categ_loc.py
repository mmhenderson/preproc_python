import numpy as np
import os, sys
import subprocess
import pandas as pd
import scipy
import scipy.io


# just make sure we're loading the right version of these modules
codepath = '/home/lab/hendersonlab/code_featsynth/preproc_code/'
# codepath = '/lab_data/hendersonlab/preproc_code/'
sys.path.insert(0, codepath)

from preproc_python import file_utils

# project_root = '/user_data/mmhender/data_UW/'
project_root = '/home/lab/hendersonlab/data_featsynth/'
retino_path = '/home/lab/hendersonlab/retino_data/ANAT/'

run_type='fLoc'


def analyze_loc(subject, subject_FS):
    
    get_task_timing(subject)
    
    run_glm_firstlevel(subject)
    
    run_glm_higherlevel(subject)
    
    organize_outputs(subject)

    make_surfaces(subject, subject_FS)
    

def get_task_timing(subject):
    
    preproc_folder = os.path.join(project_root, 'DataPreproc', subject)

    all_info_fn = os.path.join(preproc_folder, 'run_info_allsess.csv')
    run_info_allsess = pd.read_csv(all_info_fn)

    # Now loading the actual behavioral (.mat) data files in and figuring out timing.

    for ss in np.unique(np.array(run_info_allsess['session'])):

        # first i'm going to count how many nifti files (functional mri data files)
        # I expect to have for this task, for each session.
        type_names = np.array(run_info_allsess['run_type'])
        inds_this_task = np.array([run_type in t for t in type_names])
        inds_this_sess = (np.array(run_info_allsess['session'])==ss) & inds_this_task

        # these are indices of runs within the session.
        # Remember that they are one-indexed, not zero-indexed.
        run_inds = np.array(run_info_allsess['num_in_session'].iloc[inds_this_sess])
        print('\nSession %d: Functional data %s files for runs:'%(ss, run_type))
        print(run_inds)
        print('nTRs:')
        print(run_info_allsess['nTRs'].iloc[inds_this_sess])

        behav_folder = os.path.join(project_root, 'DataBehavior',subject,'Sess%01d'%ss)

        files = os.listdir(behav_folder)

        # we are looking for folder names that have "run_type" in them
        # check that this works, may need to edit depending on how files are named.
        # for fLoc: files with "_new" are the ones created by my code "convert_floc_data.m"
        # these should be load-able in python, whereas original ones aren't.
        files = [f for f in files if (run_type in f and '_new' in f)]

        if len(files)==0:
            # check if it's in a subfolder, as opposed to top-level folder
            # floc code spits out the files into subfolders by default
            subfolder = os.listdir(behav_folder)
            subfolder = [s for s in subfolder if not os.path.isfile(os.path.join(behav_folder, s))]
            for s in subfolder:
                files = os.listdir(os.path.join(behav_folder, s))
                files = [f for f in files if (run_type in f and '_new' in f)]
                files = [os.path.join(s, f) for f in files]
                if len(files)>0:
                    break

        print(files, len(files))

        # if (subject=='S01B') and (ss==2):
            
        #     # this subject we accidentally did multiple tasks, so there are 2 files.
        #     # Going to load both and combine them
        #     twoback = [f for f in files if '_2back' in f][0]
        #     oneback = [f for f in files if '_1back' in f][0]
            
        #     fn2load = os.path.join(behav_folder, oneback)
        #     print('\nLoading from %s'%fn2load)
        #     session = file_utils.load_mat_behav_data(fn2load,'session')[0]
        #     print('Found %s runs in this file'%(session['runs_done']))

        #     fn2load = os.path.join(behav_folder, twoback)
        #     print('\nLoading from %s'%fn2load)
        #     session_twoback = file_utils.load_mat_behav_data(fn2load,'session')[0]
        #     print('Found %s runs in this file'%(session_twoback['runs_done']))
        #     twoback_runsdone = session_twoback['runs_done']
        #     # ^ only want info about the runs that we did

        #     # i am making a new array that combines runs from each.
        #     session['sequence']['stim_names']=np.concatenate([np.array(session['sequence']['stim_names'])[:,None], \
        #                                                        np.array(session_twoback['sequence']['stim_names'])[:,0:twoback_runsdone]], \
        #                                                       axis=1)
        #     session['sequence']['stim_onsets']=np.concatenate([np.array(session['sequence']['stim_onsets'])[:,None], \
        #                                                        np.array(session_twoback['sequence']['stim_onsets'])[:,0:twoback_runsdone]], \
        #                                                       axis=1)

        #     n_runs = session['runs_done'] + session_twoback['runs_done']
       
        # else:

        if len(files)==0:
            print('Zero runs of %s found for Sess %d, continuing...'%(run_type, ss))
            continue
            
        assert(len(files)==1)
        fn2load = os.path.join(behav_folder, files[0])

        print('\nLoading from %s'%fn2load)
        session = file_utils.load_mat_behav_data(fn2load,'session')[0]

        # this marks how many runs actually completed. even if you made "sequence" for more than this number.
        n_runs = session['runs_done']

              
        # should be same number as we have brain data for.
        assert(n_runs==len(run_inds))

            
        print('Session %d: Found %d total runs in behavior files.'%(ss, n_runs))

        # loop over my runs, within this session
        for ri, rr in enumerate(run_inds):

            # Getting a list of all the images that were shown, and their onsets in sec
            # these are individual images within each block.
            if len(np.array(session['sequence']['stim_names']).shape)>1:
                # if >1 run, this is an array [n_images x n_runs]
                image_list = np.array(session['sequence']['stim_names'])[:,ri]
                stim_onsets = np.array(session['sequence']['stim_onsets'])[:,ri]
            else:
                # else it is just one list
                image_list = np.array(session['sequence']['stim_names'])
                stim_onsets = np.array(session['sequence']['stim_onsets'])

            # adjust for time of the countdown period
            count_down = session['count_down']
            stim_onsets = stim_onsets + count_down

            # Then figure out condition of each image.
            # Assumes your image names are formatted like <cond>-<number>
            cond_list = np.array([i.split('-')[0] for i in image_list])

            # now figure out where condition is changing. 
            cond_names, cond_inds = np.unique(cond_list, return_inverse=True)
            print('\nRun %d: Conditions are:'%rr)
            print(cond_names)
            start_block_inds = (np.diff(np.concatenate([np.array([100]), cond_inds])))!=0
            # NOTE that sometimes, due to randomization of image sequence, 
            # we will have two blocks in a row of same condition. This is on purpose and ok.
            # For the purpose of our analysis, will treat these as one big long block.


            block_onsets = stim_onsets[start_block_inds]

            block_cond_inds = np.array(cond_inds)[start_block_inds]
            block_cond_names = np.array(cond_names)[block_cond_inds]

            # when stimuli stopped being presented, end of the whole block
            last_block_end = session['total_run_dur']

            print('Run is %d seconds (double check this makes sense with your TRs!)'%session['total_run_dur'])

            # how long is each block? again this includes any double-length.
            block_durs = np.diff(np.concatenate([block_onsets, np.array([last_block_end])]))

            # here I am double checking the timing. Blocks should all be expected length or 2x that length.
            block_len_expected = session['sequence']['stim_per_block'] * session['sequence']['stim_duty_cycle']
            assert(np.all(np.mod(block_durs, block_len_expected)==0))

            # Now write the text files.
            output_dir = os.path.join(project_root, 'CategLocAnalysis', subject, 'EVs');
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)

            # Going to make one text file per condition.
            # Where conditions are like: face, place, baseline, etc.
            # They are formatted as [onset, duration, condition (1/0 binary)]
            # We call them "EV" files, and they will get fed into FSL's FEAT program.
            for cond_name in cond_names:

                text_file = os.path.join(output_dir, '%s_Sess%02d_Run%02d_%s.txt'%(subject, ss, rr, cond_name))
                print('Writing to %s'%text_file)

                with open(text_file,'w') as h:

                    # make 0/1 for this condition or not
                    this_cond_int = (block_cond_names==cond_name).astype(int)
                    for xx in range(len(block_onsets)):
                        line = '%.1f  %.1f  %d\n'%(block_onsets[xx], block_durs[xx], this_cond_int[xx])
                        h.write(line)

                        
                        
def run_glm_firstlevel(subject, spatial_smoothing = False):

    preproc_folder = os.path.join(project_root, 'DataPreproc', subject)

    all_info_fn = os.path.join(preproc_folder, 'run_info_allsess.csv')
    run_info_allsess = pd.read_csv(all_info_fn)

    # these are the 8 conditions in my experiment.
    # I am hard-coding these so that the order is the same
    # but will need to change for future experiments.
    cond_names = ['adult', 'baseline', 'body', 'car', \
                  'food', 'house', 'instrument','scrambled']

    n_preds = 8

    ev_dir = os.path.join(project_root, 'CategLocAnalysis', subject, 'EVs');
    fsf_dir = os.path.join(project_root, 'CategLocAnalysis', subject, 'fsfs');
    feat_dir = os.path.join(project_root, 'CategLocAnalysis', subject, 'feats');
    if not os.path.exists(fsf_dir):
        os.makedirs(fsf_dir)
    if not os.path.exists(feat_dir):
        os.makedirs(feat_dir)

    nifti_dir = os.path.join(project_root, 'DataPreproc',subject)
    all_info_fn = os.path.join(preproc_folder, 'run_info_allsess.csv')
    run_info_allsess = pd.read_csv(all_info_fn)

    # Start from this template, made using the FEAT dialog box.
    # Might need to edit for a new experiment with diff number of conditions.
    template_file = os.path.join(project_root, 'CategLocAnalysis', 'templates', 'design_template.fsf')

    # find all runs in this task
    type_names = np.array(run_info_allsess['run_type'])
    inds_this_task = np.where(np.array([run_type in t for t in type_names]))[0]

    # Looping over individual runs.
    # This is the "lower-level" FEAT analysis.
    for ii in inds_this_task:

        # this is the most fully pre-processed file, use as input to FEAT.
        nifti_fn = os.path.join(preproc_folder, run_info_allsess['reg_mc_det_fn'].iloc[ii])    

        ss = run_info_allsess['session'].iloc[ii]
        rr = run_info_allsess['num_in_session'].iloc[ii]
        n_trs = run_info_allsess['nTRs'][ii]
        tr_length = run_info_allsess['TR_length'][ii]

        print('For session %d run %d, we have %d TRs of length %.2fsec. Check this is what you expect!'%(ss, rr, n_trs, tr_length))

        # define where files will go for this run
        output_dir = os.path.join(feat_dir, '%s_Sess%02d_Run%02d'%(subject, ss, rr))

        output_dir_actual = output_dir+'.feat'
        if os.path.exists(output_dir_actual):
            
            print('Directory already exists at: %s'%(output_dir_actual))
            print('NOT running analysis again since directory exists. Delete this directory if you want to re-do it.')
            print('Skipping to next run...\n')
            continue
            
        else:

            with open(template_file, 'r') as f:

                lines_template = f.readlines()

            # Now i'm going to run through the lines in this text file, and change things.
            # This allows us to update the template based on individual runs
            # I'm not making any huge changes here like a change in number of conditions 
            # for example. To make bigger changes, make a new template using FEAT dialog box.
            lines_new = []
            for l in lines_template:

                lnew = l

                # directory where all results get saved
                if 'set fmri(outputdir)' in l:
                    lnew = 'set fmri(outputdir) "%s"\n'%(output_dir)
                    print(lnew)

                # input nifti functional data file, for this run
                if 'set feat_files(1)' in l:
                    lnew = 'set feat_files(1) "%s"\n'%(nifti_fn)
                    print(lnew)

                # set number of TRs
                if 'set fmri(npts)' in l:
                    lnew = 'set fmri(npts) %d\n'%(n_trs)
                    print(lnew)

                if 'set fmri(tr)' in l:
                    lnew = 'set fmri(tr) %.4f\n'%(tr_length)
                    print(lnew)

                # spatial smoothing on or off
                if 'set fmri(smooth)' in l:
                    if spatial_smoothing:
                        lnew = 'set fmri(smooth) 2.0\n'
                    else:
                        lnew = 'set fmri(smooth) 0\n'
                    print(lnew)

                # files that define the time sequence of events in the run 
                # (EV files, which we made earlier)
                for cc, cond_name in zip(np.arange(1, len(cond_names)+1), cond_names):

                    ev_name = os.path.join(ev_dir, '%s_Sess%02d_Run%02d_%s.txt'%(subject, ss, rr, cond_name))

                    if 'set fmri(custom%d)'%(cc) in l:
                        lnew = 'set fmri(custom%d) "%s"\n'%(cc, ev_name)
                        print(lnew)

                    if 'set fmri(evtitle%d)'%(cc) in l:
                        print(l)
                        # double check that the cond names are correct here
                        # this is important because contrasts are set based on these numbers, 
                        # and we want to make sure they're correct.
                        assert(cond_name in l)

                lines_new += [lnew]


            # put the edited text into new file.
            # This becomes the input to FSL FEAT
            new_design_file = os.path.join(fsf_dir, '%s_Sess%02d_Run%02d.fsf'%(subject, ss, rr))
            print('Writing design to %s\n\n'%new_design_file)

            with open(new_design_file, 'w') as f:
                f.writelines(lines_new)



            print('Starting GLM (FEAT analysis) for Session %d, Run %d\n (this may take a while...)'%(ss, rr))
            # now we're using the design file as input, and running the analysis
            cmd = 'feat %s'%(new_design_file)
            print(cmd)
            sys.stdout.flush()
            err = subprocess.call(cmd, shell=True)
            if err:
                print('Error in FEAT analysis for Session %d, Run %d'%(ss, rr))
            else:
                print('Done with Session %d, Run %d'%(ss, rr))

            # Make some extra matrices as a workaround that enables higher-level analysis. 
            # The FEAT program expects there to be a registration of your data between native space
            # and standardized (like MNI) space). We don't need this because we're doing this in native 
            # space. So I'm putting an identity matrix in these files. 
            # https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FEAT/UserGuide

            reg_dir = os.path.join(output_dir_actual, 'reg')
            if not os.path.exists(reg_dir):
                os.makedirs(reg_dir)

            reg_placeholder = np.identity(4).astype(int)

            fn = os.path.join(reg_dir, 'example_func2standard.mat')
            print('\nWriting to %s'%fn)
            file_utils.write_xfm_mat(reg_placeholder, fn)

            fn = os.path.join(reg_dir, 'standard2example_func.mat')
            print('Writing to %s'%fn)
            file_utils.write_xfm_mat(reg_placeholder, fn)

            # copy example here as my "standard" template
            example_file = os.path.join(output_dir_actual, 'example_func.nii.gz')
            cmd = 'scp %s %s'%(example_file, os.path.join(reg_dir, 'standard.nii.gz'))
            print(cmd)
            sys.stdout.flush()
            err = subprocess.call(cmd, shell=True)

        
def run_glm_higherlevel(subject, higher_level_FE = True):

    # Now prepare for the "higher-level" analysis, which combines results across runs.

    preproc_folder = os.path.join(project_root, 'DataPreproc', subject)

    all_info_fn = os.path.join(preproc_folder, 'run_info_allsess.csv')
    run_info_allsess = pd.read_csv(all_info_fn)

    # find all runs in this task
    type_names = np.array(run_info_allsess['run_type'])
    inds_this_task = np.where(np.array([run_type in t for t in type_names]))[0]
    # ^ inds_this_task are indices into the big df with all runs
    n_runs = len(inds_this_task)
    
    tr_length = run_info_allsess['TR_length'][inds_this_task[0]]
   
    run_inds_simple = np.arange(1, n_runs+1)
    # ^ this is a simple index of runs in this task only, start from 1
    
    # This is the template we're going to edit.
    template_file = os.path.join(project_root, 'CategLocAnalysis', 'templates', 'high_level_design_template_%d.fsf'%(n_runs))

    # define where files will go for this analysis
    fsf_dir = os.path.join(project_root, 'CategLocAnalysis', subject, 'fsfs');
    feat_dir = os.path.join(project_root, 'CategLocAnalysis', subject, 'feats');
    
    output_dir = os.path.join(feat_dir, 'AllSessionsFE')

    # load the template and get ready to edit.
    with open(template_file, 'r') as f:
        lines = f.readlines()

    lines_new = []
    for l in lines:

        lnew = l

        # Here, npts is number of input runs. NOT TRs
        if 'set fmri(npts)' in l:
            print(l)
            n_runs_expected = int(l.split('set fmri(npts) ')[1].split('\n')[0])
            assert(n_runs_expected==n_runs)

        if 'set fmri(tr)' in l:
            lnew = 'set fmri(tr) %.4f\n'%(tr_length)
            print(lnew)

        # which method will we use to combine across runs?
        if 'set fmri(mixed_yn)' in l:
            print(l)
             # '# Higher-level modelling\n',
             # '# 3 : Fixed effects\n',
             # '# 0 : Mixed Effects: Simple OLS\n',
             # '# 2 : Mixed Effects: FLAME 1\n',
             # '# 1 : Mixed Effects: FLAME 1+2\n',
            if higher_level_FE:
                lnew = 'set fmri(mixed_yn) 3\n'
            else:
                lnew = 'set fmri(mixed_yn) 2\n'
            print(lnew)

        if 'set fmri(outputdir)' in l:
            lnew = 'set fmri(outputdir) %s\n'%(output_dir)
            print(lnew)

        # adding input data names for each run
        for ri, ii in zip(run_inds_simple, inds_this_task):

            # print(ri, ii, ss, rr)
            ss = run_info_allsess['session'].iloc[ii]
            rr = run_info_allsess['num_in_session'].iloc[ii]

            if 'set feat_files(%d)'%ri in l:
                print(l)
                run_output_dir_actual = os.path.join(feat_dir, '%s_Sess%02d_Run%02d.feat'%(subject, ss, rr))
                lnew = 'set feat_files(%d) "%s"\n'%(ri, run_output_dir_actual)
                print(lnew)


        lines_new += [lnew]

    # Now writing this to a new file.
    # This becomes the input to FSL FEAT
    new_design_file = os.path.join(fsf_dir, '%s_high_level_design.fsf'%subject)
    print('Writing design to %s\n\n'%new_design_file)

    with open(new_design_file, 'w') as f:
        f.writelines(lines_new)

    # Here is where we will run the feat program itself
    output_dir_actual = output_dir+'.gfeat'

    if os.path.exists(output_dir_actual):
        
        print('Directory already exists at: %s'%(output_dir_actual))
        print('NOT running analysis again since directory exists. Delete this directory if you want to re-do it.')
        
    else:

        print('\nRunning GLM (FEAT analysis) across all runs')

        cmd = 'feat %s'%(new_design_file)
        print(cmd)
        sys.stdout.flush()
        err = subprocess.call(cmd, shell=True)
        if err:
            print('\nError in FEAT analysis when combining across runs')

        else:
            print('\nDONE')
            
            
def organize_outputs(subject):
    
    # this is just a convenience function to help put all our FEAT outputs into one folder
    # we're just copying and renaming many files.
    # Helps avoid a lot of clicking when you're loading these files into FSLeyes 

    
    feat_dir = os.path.join(project_root, 'CategLocAnalysis', subject, 'feats');

    # this is where i'm putting things
    
    outfolder = os.path.join(project_root, 'CategLocAnalysis', subject, 'summary_stats')
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)

    print('\nCopying files to %s...'%outfolder)
    
    fdirs = os.listdir(feat_dir)
    fdirs = [d.split('.')[0] for d in fdirs]

    fdirs_avg = [d for d in fdirs if 'AllSessions' in d]
    fdirs_ss = [d for d in fdirs if 'AllSessions' not in d]
   
    # load first FSF file to figure out contrast names, etc
    fsf_dir = os.path.join(project_root, 'CategLocAnalysis', subject, 'fsfs');

    fn = os.path.join(fsf_dir, '%s.fsf'%(fdirs_ss[0]))

    with open(fn) as f:
        lines = f.readlines()

    ind = [li for li, l in enumerate(lines) if l=='# Number of contrasts\n']
    l = lines[ind[0]+1]
    n_contrasts = int(l.split('set fmri(ncon_orig) ')[1].split('\n')[0])
  
    ind = np.array([li for li, l in enumerate(lines) if 'Title for contrast_real' in l])
    lines_use = np.array(lines)[ind+1]
    cnames = [l.split('"')[1] for l in lines_use]
    cnames = [c.replace('>', '_gr_') for c in cnames]
   
    assert(len(cnames)==n_contrasts)
    
    for si, fdir in enumerate(fdirs_ss):
        statsfolder = os.path.join(feat_dir, '%s.feat'%fdir, 'stats')
        for cc, cname in zip(np.arange(1, n_contrasts+1), cnames):
            statsfilename = os.path.join(statsfolder, 'tstat%d.nii.gz'%(cc))
            newname = os.path.join(outfolder, '%s_%s_tstat.nii.gz'%(fdir, cname))
            cmd = 'cp %s %s'%(statsfilename, newname)
            print(cmd)
            sys.stdout.flush()
            err = subprocess.call(cmd, shell=True)


    for si, fdir in enumerate(fdirs_avg):

        for cc, cname in zip(np.arange(1, n_contrasts+1), cnames):

            statsfolder = os.path.join(feat_dir, '%s.gfeat'%fdir, 'cope%d.feat'%cc,'stats')

            statsfilename = os.path.join(statsfolder, 'tstat1.nii.gz')
            newname = os.path.join(outfolder, '%s_%s_%s_tstat.nii.gz'%(subject, fdir, cname))
            cmd = 'cp %s %s'%(statsfilename, newname)
            print(cmd)
            sys.stdout.flush()
            err = subprocess.call(cmd, shell=True)



def make_surfaces(subject, subject_FS, debug=False):

    # Projecting my t-statistic maps, which are in 3D volume space originally, into surface space
    # These files can then be loaded in freeview to plot the data on inflated surface.

    
    preproc_folder = os.path.join(project_root, 'DataPreproc', subject)
    
    all_info_fn = os.path.join(preproc_folder, 'run_info_allsess.csv')
    
    # this is where outputs of GLM are saved (t-stats)
    input_folder = os.path.join(project_root, 'CategLocAnalysis', subject, 'summary_stats')
    inputs = os.listdir(input_folder)
    # inputs = [i for i in inputs if (('AllSessions' in i) and ('proj' not in i))]
    inputs = [i for i in inputs if (('proj' not in i))]
    
    surf = 'white'
    # surf = 'pial'
    
    # loop over these inputs, and make surface files that we can load into freeview
    for i in inputs:
    
        input_vol_file = os.path.join(input_folder, i)
    
        # do this for each hemisphere separately
        for hemi in ['lh','rh']:            
            
            new = i.split('.nii.gz')[0] + '_proj_%s_%s'%(hemi, surf) + '.nii'
            
            proj_surf_file = os.path.join(input_folder, new)
    
            if os.path.exists(proj_surf_file):
                print('%s already exists, skipping'%(proj_surf_file))
                continue
    
            # Need this file: it tells us how to map from functional to anatomical space
            # Made during preprocessing.py
            reg_file = os.path.join(preproc_folder, 'Sess01_func2anat_xfm.dat')
    
            # project to surface 
            cmd = 'mri_vol2surf ' + \
                  '--src %s '%(input_vol_file) + \
                  '--out %s '%(proj_surf_file) + \
                  '--srcreg %s '%(reg_file) + \
                  '--out_type .surf ' + \
                  '--hemi %s --surf %s'%(hemi, surf)
            
            print(cmd)
            sys.stdout.flush()
            if not debug:
                err = subprocess.call(cmd, shell=True)
