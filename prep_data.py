import numpy as np
import os, sys
import subprocess
import pandas as pd
from datetime import datetime
import nibabel as nib


def organize_functional(inpath):

    # This function will take your raw data and help set up the file structure that we want.
    # all data should start out in your "Dicoms" folder, and we're going to move some to a 
    # "Niftis" folder so we can easily find it later. 
    # If the niftis don't exist yet, we'll also try to create them here.
    
    # inpath should be something like: /user_data/mmhender/data/DataRaw/S01/Session1/
    # which has "Dicoms" and "Niftis" inside it.
    
    dicom_dir = os.path.join(inpath, 'Dicoms')
    nifti_dir = os.path.join(inpath, 'Niftis')
    
    if not os.path.exists(nifti_dir):
        os.makedirs(nifti_dir)
    
    # list all folders - this is functional, anatomical, everything.
    raw_folders = os.listdir(dicom_dir)
    
    # just grabbing functional runs and field maps here
    raw_folders = [f for f in raw_folders if \
                   (('func' in f or 'fmap' in f) and 'Physio' not in f)]
    
    for f in raw_folders:
    
        print(f)
    
        done=False
        nii_dir_expected = os.path.join(nifti_dir, f)
        if os.path.exists(nii_dir_expected):
            files = os.listdir(nii_dir_expected)
            nii_file = [f for f in files if '.nii.gz' in f]
            if len(nii_file)>0:
                done=True
                print('done')
    
        if not done:
     
            # seeing what files are already here
            this_dir = os.path.join(dicom_dir, f)
            files = os.listdir(this_dir)
    
            nii_file = [f for f in files if '.nii.gz' in f]
            print(nii_file)
        
            if len(nii_file)==0:
                
                # in this case we need to create the nifti file from the dicom files
                # usually for flywheel data this is already done.
                
                # dfolder is the folder of dicom images, this would be the input to dcm2niix
                dfolder = [f for f in files if ((f[-6:]=='.dicom') and \
                           os.path.isdir(os.path.join(this_dir, f)))]
                print(dfolder)
                # if this folder doesn't exist yet, we will need to unzip the zipped files
                if len(dfolder)==0:    
                    dfile_zip = [f for f in files if 'dicom.zip' in f]
                    for df in dfile_zip:
                        fn2unzip = os.path.join(this_dir, df)
                        cmd = 'unzip %s -d %s'%(fn2unzip, this_dir)
                        print(cmd)
                        err = subprocess.call(cmd, shell=True)
                        if err:
                            print('Error in dicom unzipping')
        
                    # get path to the new unzipped dir
                    files = os.listdir(this_dir)
                    dfolder = [f for f in files if ((f[-6:]=='.dicom') and \
                               os.path.isdir(os.path.join(this_dir, f)))]
                   
                dicom_folder = os.path.join(this_dir, dfolder[0])
    
                # now do the actual conversion.
                # make sure this function dcm2niix is on your path!
                cmd = 'dcm2niix -o %s -z y -f '%this_dir + '%' + 'f %s'%dicom_folder
                print(cmd)
                err = subprocess.call(cmd, shell=True)
    
            # Now we have the nifti file. 
            # Going to move it now just so it better aligns with desired file structure.
            files = os.listdir(this_dir)
            nii_file = [f for f in files if ('.nii.gz' in f)]
            print(nii_file)
            
            new_dir = os.path.join(nifti_dir, f)
            if not os.path.exists(new_dir):
                os.makedirs(new_dir)
            for nf in nii_file:
                # only move the file if it's not moved yet 
                if not os.path.exists(os.path.join(new_dir, nf)):
                    cmd = 'mv %s %s'%(os.path.join(this_dir, nf), \
                                      os.path.join(new_dir, nf))
                    print(cmd)
                    err = subprocess.call(cmd, shell=True)
                    if err:
                        print('Error in moving niftis')
                else:
                    print('skipping move for this file')
                    
                    
def get_run_info(inpath):
    
    # This is a function that goes into your "Niftis" folder and makes a pd.DataFrame
    # Saves it as a .csv file that has key info about the runs.
    # Once this is done, you should look over the .csv file and make sure everything aligns 
    # with what you are expecting. 
    # if there's anything weird, move the bad niftis manually to a new sub-folder, then 
    # run this again to remake the .csv file.

    nifti_dirs = os.listdir(os.path.join(inpath, 'Niftis'))
    nifti_dirs = [n for n in nifti_dirs if ('func' in n or 'fmap' in n)]
    run_info = pd.DataFrame(columns=['num_in_session','dir_name', 'run_type','run_number',\
                                     'date_raw','dstring',\
                                     'time_raw', 'tstring', \
                                     'size_mb', \
                                     'nifti_fn', \
                                     'nTRs', 'TR_length'])

    for di, nd in enumerate(nifti_dirs):

        thisdir = os.path.join(inpath, 'Niftis', nd)

        # gathering info from the file name.
        # these are somewhat specific to my experiment here, might want to check 
        # for other new experiments.
        if 'func' in nd:
            run_type = nd.split('func-bold_task-')[1].split('_run')[0]
            if 'run' in nd:
                run_number = int(nd.split('_run-')[1][0:2])
        else:
            if 'AP' in nd:
                run_type = 'fmap_AP'
            elif 'PA' in nd:
                run_type = 'fmap_PA'
            run_number = np.nan


        fnames = os.listdir(thisdir)
        fnames = [f for f in fnames if (f[0]=='1' and '.nii.gz' in f and 'topup' not in f)]
        print(thisdir)
        print(fnames)
        assert(len(fnames)==1)
        n = fnames[0]

        full_name = os.path.join(thisdir, n)
        # print(full_name)

        # Get approx size. This helps us know if something is really wrong.
        # For example if a scan got stopped in the middle, it will be small.
        # Or if it is corrupted in some way, it will also be small.
        size_mb = np.round(os.path.getsize(full_name)/(1024**2), 1)
        
        if not np.isnan(run_number):
            nfile = nib.load(full_name)
            dat = np.array(nfile.get_fdata())
            # Number of TRs, last dimension here
            nTRs = dat.shape[3]
            # Find TR length. 
            # Pixdim elements [1,2,3,4] are [x,y,z,t]
            tr_length = nfile.header['pixdim'][4]
            # CHECK both these against what you're expecting it to be for this protocol.
        else:
            nTRs = np.nan
            tr_length = np.nan
            
        # Figure out run acquisition time, based on string in the file name.
        traw = n.split('_')[9]

        date_raw = traw[0:8]
        dstring = datetime.strftime(datetime.strptime(date_raw, '%Y%m%d'), '%D')

        time_raw = traw[8:14]
        tstring = datetime.strftime(datetime.strptime(time_raw, '%H%M%S'), '%I:%M %p')

        # Putting these elements into our big dataframe, one run at a time
        run_info = run_info.append(pd.DataFrame({'dir_name': nd, \
                                                 'run_type': run_type, \
                                                 'run_number': run_number, \
                                                 'date_raw': date_raw,
                                                 'dstring': dstring, \
                                                 'time_raw': time_raw,
                                                 'tstring': tstring, \
                                                 'size_mb': size_mb, \
                                                 'nifti_fn': full_name, \
                                                 'nTRs': nTRs, \
                                                 'TR_length': tr_length, 
                                                }, \
                                               index = [di]))

    # now sorting by time within this session.

    sorder = np.argsort(run_info['time_raw'])
    run_info = run_info.iloc[sorder]

    run_info['num_in_session'] = np.arange(1, run_info.shape[0]+1)

    run_info = run_info.set_index(np.arange(run_info.shape[0]))

    save_path = os.path.join(inpath, 'Niftis', 'run_info.csv')
    print('saving to %s'%save_path)
    run_info.to_csv(save_path, index=False)


def process_tar_raw(inpath):

    # this function will take the raw .tar file that you download from flywheel
    # (using fw download [path] with the CLI, or something similar)
    # we will unzip the .tar file, find where all the folders for individual 
    # scans are located, and place all those folders in a folder called "Dicoms"
    # which goes inside your "inpath"

    tarfile = os.listdir(inpath)
    tarfile = [f for f in tarfile if '.tar' in f]
    assert(len(tarfile)==1)
    tarfile_full = os.path.join(inpath, tarfile[0])

    tmpname = tarfile[0].split('.tar')[0] + '_tmp'
    tmpfolder = os.path.join(inpath, tmpname)
    if not os.path.exists(tmpfolder):
        os.makedirs(tmpfolder)

    # Unzip the .tar file here
    print('Unzipping data from:\n%s to\n%s'%(tarfile_full, tmpfolder))
    cmd = 'tar -xvf %s -C %s'%(tarfile_full, tmpfolder)
    print(cmd)
    err = subprocess.call(cmd, shell=True)
    
    # figure out where my data lives in here
    path_to_dat = tmpfolder
    for ii in range(10):
        f = os.listdir(path_to_dat)
        if len(f)==1:
            path_to_dat = os.path.join(path_to_dat, f[0])
        else:
            # this is where my subdirectories live
            continue
    
    print('Data is in: %s'%path_to_dat)
    
    print('Files are:')
    print(os.listdir(path_to_dat))
    
    # Now i'm moving all the subdirectories into a folder called "Dicoms"
    new_path = os.path.join(inpath, 'Dicoms/')
    if not os.path.exists(new_path):
        os.makedirs(new_path)
    
    print('Moving all folders to %s'%(new_path))
    
    cmd = 'mv %s/* %s'%(path_to_dat, new_path)
    print(cmd)
    err = subprocess.call(cmd, shell=True)
 
    
    # finally, fix any folders that have a space in the name
    # this can create bugs later on.
    for subfolder in os.listdir(new_path):
        if ' ' in subfolder:
            new_subfolder = subfolder.replace(' ', '_')
            old_subfolder = subfolder.replace(' ', '\ ')
            cmd = 'mv %s %s'%(os.path.join(new_path, old_subfolder), \
                            os.path.join(new_path, new_subfolder))
            print(cmd)
            err = subprocess.call(cmd, shell=True)

   
    if err==0:
        cmd = 'rm -r %s'%(tmpfolder)
        print(cmd)
        err = subprocess.call(cmd, shell=True)
