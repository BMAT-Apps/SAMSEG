#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 14:04:17 2024

@author: colin
"""

import os
from os.path import join as pjoin
from os.path import exists as pexists
import sys
import argparse
import time
import shutil
import subprocess
import nibabel as nib


def bids_SAMSEG(bids, sub, ses, mprage=None, flair=None, cpus=8):    
    
    if mprage == None or mprage == "":
        mprage_bool = False
    else:
        mprage_bool = True
    
    if flair == None or flair == "":
        flair_bool = False
    else:
        flair_bool = True
    
    segment = 'segmentation'
    transfo = 'transformation'
    sub_ses_anat = pjoin(bids, f'sub-{sub}', f'ses-{ses}', 'anat')
    flair = f'sub-{sub}_ses-{ses}_{flair}.nii.gz'
    mprage = f'sub-{sub}_ses-{ses}_{mprage}.nii.gz'
    sub_ses_segment = pjoin(bids, 'derivatives', 'SAMSEG', f'sub-{sub}', f'ses-{ses}', segment)
    sub_ses_transfo = pjoin(bids, 'derivatives', 'SAMSEG', f'sub-{sub}', f'ses-{ses}', transfo)
    
    # Perform Lesion Segmentation, Localisation and Volumetry        
    print('Start Preprocessing...')
    
    # Pre check
    ## Create corresponding directory if necessary
    os.makedirs(sub_ses_segment, exist_ok=True)
    os.makedirs(sub_ses_transfo, exist_ok=True)
    
    # Pre check
    ## check if preprocessing has already been done previously, and delete the results if necessary
    if os.path.isdir(pjoin(sub_ses_transfo, f'sub-{sub}_ses-{ses}_FLAIR_preprocessed')):
        os.remove(pjoin(sub_ses_transfo, f'sub-{sub}_ses-{ses}_FLAIR_preprocessed'))
    if pexists(pjoin(sub_ses_transfo, f'sub-{sub}_ses-{ses}_MPRAGE_preprocessed.nii.gz')):
        os.remove(pjoin(sub_ses_transfo, f'sub-{sub}_ses-{ses}_MPRAGE_preprocessed.nii.gz'))
          
    if not pexists(pjoin(sub_ses_anat, flair)):
        flair_bool = False
        print(f'No FLAIR for sub-{sub} ses-{ses}')
    if not pexists(pjoin(sub_ses_anat, mprage)):
        mprage_bool = False
        print(f'No MPRAGE for sub-{sub} ses-{ses}')
        
    if mprage_bool and flair_bool:
        shutil.copyfile(pjoin(sub_ses_anat, mprage), pjoin(sub_ses_transfo, f'sub-{sub}_ses-{ses}_MPRAGE.nii.gz'))
        
        try:
            print(f'Registration of FLAIR on MPRAGE for sub-{sub} ses-{ses}...')
            subprocess.Popen(f'mri_coreg --mov {pjoin(sub_ses_anat, flair)} --ref {pjoin(sub_ses_anat, mprage)} --reg {pjoin(sub_ses_transfo, "flairToT1.lta")}', shell=True).wait()
            subprocess.Popen(f'mri_vol2vol --mov {pjoin(sub_ses_anat, flair)} --reg {pjoin(sub_ses_transfo, "flairToT1.lta")} --o {pjoin(sub_ses_transfo, f"sub-{sub}_ses-{ses}_FLAIR.nii.gz")} --targ {pjoin(sub_ses_anat, mprage)}', shell=True).wait()                                
            
        except Exception as e:
            print(f'[ERROR] - SAMSEG | {e} during registration of FLAIR on MPRAGE for sub-{sub}_ses{ses}!')
            return
        
    elif flair_bool:
        shutil.copyfile(pjoin(sub_ses_anat, flair), pjoin(sub_ses_transfo, f'sub-{sub}_ses-{ses}_FLAIR.nii.gz'))
    elif mprage_bool:
        shutil.copyfile(pjoin(sub_ses_anat, mprage), pjoin(sub_ses_transfo, f'sub-{sub}_ses-{ses}_MPRAGE.nii.gz'))
    else:
        return
            
        print('End Preprocessing!')
    
    ## Step2 : Segmentation using SAMSEG to compute lesion probability mask
        
    print('Start Segmentating lesions with SAMSEG...')
    
    # Actions
    ## Run SAMSEG with FLAIR and MPRAGE
    if mprage_bool and flair_bool:
        
        print('Segmenting with FLAIR & mprage')
        
        # Actions
        try:
            print(f'Running Samseg for sub-{sub} ses-{ses}...')
            
            subprocess.Popen(f"mri_convert {pjoin(sub_ses_transfo, f'sub-{sub}_ses-{ses}_FLAIR.nii.gz')} {pjoin(sub_ses_transfo, f'sub-{sub}_ses-{ses}_FLAIR.mgz')}", shell=True).wait()
            
            subprocess.Popen(f"mri_convert {pjoin(sub_ses_transfo, f'sub-{sub}_ses-{ses}_MPRAGE.nii.gz')} {pjoin(sub_ses_transfo, f'sub-{sub}_ses-{ses}_MPRAGE.mgz')}", shell=True).wait()
            
            subprocess.Popen(f"run_samseg -i {pjoin(sub_ses_transfo, f'sub-{sub}_ses-{ses}_MPRAGE.mgz')} -i {pjoin(sub_ses_transfo, f'sub-{sub}_ses-{ses}_FLAIR.mgz')} --lesion --lesion-mask-pattern 0 1 --threads {cpus} -o {pjoin(sub_ses_segment,'SAMSEG_results')} --save-posteriors", shell=True).wait()
            
            subprocess.Popen(f"mri_convert {pjoin(sub_ses_segment, 'SAMSEG_results', 'posteriors', 'Lesions.mgz')} {pjoin(sub_ses_segment, f'sub-{sub}_ses-{ses}_lesions.nii.gz')}", shell=True).wait()
            
        except Exception as e:
            print(f'[ERROR] - SAMSEG | {e} while running Samseg for sub-{sub}_ses{ses}!')
            return
        
    ## Run SAMSEG with only FLAIR 
    elif flair_bool:
        
        print('Segmenting with only FLAIR')
       
        # Actions
        try:
            print(f'Running Samseg for sub-{sub} ses-{ses}...')

            subprocess.Popen(f"mri_convert {pjoin(sub_ses_transfo, f'sub-{sub}_ses-{ses}_FLAIR.nii.gz')} {pjoin(sub_ses_transfo, f'sub-{sub}_ses-{ses}_FLAIR.mgz')}", shell=True).wait()
            
            subprocess.Popen(f"run_samseg -i {pjoin(sub_ses_transfo, f'sub-{sub}_ses-{ses}_FLAIR.mgz')} --lesion --lesion-mask-pattern 0 1 --threads {cpus} -o {pjoin(sub_ses_segment,'SAMSEG_results')} --save-posteriors", shell=True).wait()
            
            subprocess.Popen(f"mri_convert {pjoin(sub_ses_segment, 'SAMSEG_results', 'posteriors', 'Lesions.mgz')} {pjoin(sub_ses_segment, f'sub-{sub}_ses-{ses}_lesions.nii.gz')}", shell=True).wait()
            
        except Exception as e:   
            print(f'[ERROR] - SAMSEG | {e} while running Samseg for sub-{sub}_ses{ses}!')
            return
        
    elif mprage_bool:
        print('Segmenting with only MPRAGE')
        
        # Actions
        try:
            print(f'Running Samseg for sub-{sub} ses-{ses}...')

            subprocess.Popen(f"mri_convert {pjoin(sub_ses_transfo, f'sub-{sub}_ses-{ses}_MPRAGE.nii.gz')} {pjoin(sub_ses_transfo, f'sub-{sub}_ses-{ses}_MPRAGE.mgz')}", shell=True).wait()
            
            subprocess.Popen(f"run_samseg -i {pjoin(sub_ses_transfo, f'sub-{sub}_ses-{ses}_MPRAGE.mgz')} --lesion --lesion-mask-pattern 0 1 --threads {cpus} -o {pjoin(sub_ses_segment,'SAMSEG_results')} --save-posteriors", shell=True).wait()
            
            subprocess.Popen(f"mri_convert {pjoin(sub_ses_segment, 'SAMSEG_results', 'posteriors', 'Lesions.mgz')} {pjoin(sub_ses_segment, f'sub-{sub}_ses-{ses}_lesions.nii.gz')}", shell=True).wait()
            
        except Exception as e:   
            print(f'[ERROR] - SAMSEG | {e} while running Samseg for sub-{sub}_ses{ses}!')
            return
        
    ## Binarizing the lesion probability  mask    
    # Pre check
    ## check if files exist
    if not pexists(pjoin(sub_ses_segment, f'sub-{sub}_ses-{ses}_lesions.nii.gz')):
        file = pjoin(sub_ses_segment, f'sub-{sub}_ses-{ses}_lesions.nii.gz')
        print(f'[ERROR] - SAMSEG | FileNotFound: File "{file}" not Found !')
        return
    
    try:                
        # Actions
        print(f'Binarizing lesion probability mask for sub-{sub} ses-{ses}...')
       
        threshold = 0.5
        image = nib.load(pjoin(sub_ses_segment, f'sub-{sub}_ses-{ses}_lesions.nii.gz'))
        lesions = image.get_fdata()
        lesions[lesions >= threshold] = 1
        lesions[lesions < threshold] = 0

        
        lesions_nifti = nib.Nifti1Image(lesions, affine=image.affine)
        
        nib.save(lesions_nifti, pjoin(sub_ses_segment, f'sub-{sub}_ses-{ses}_lesions_binary.nii.gz'))

        
    except Exception as e:
        print(f'[ERROR] - SAMSEG | {e} while binarizing lesion probability mask for sub-{sub}_ses{ses}!')
        return
    
    print('End SAMSEG!')
    
    
    
def bids_SAMSEG_docker(bids, sub, ses, mprage=None, flair=None, cpus=8):    
    
    path = os.path.dirname(os.path.abspath(__file__))
    fs_license = pjoin(path, 'license.txt')
    
    if mprage == None or mprage == "":
        mprage_bool = False
    else:
        mprage_bool = True
    
    if flair == None or flair == "":
        flair_bool = False
    else:
        flair_bool = True
    
    segment = 'segmentation'
    transfo = 'transformation'
    sub_ses_anat = pjoin(bids, f'sub-{sub}', f'ses-{ses}', 'anat')
    flair = f'sub-{sub}_ses-{ses}_{flair}.nii.gz'
    mprage = f'sub-{sub}_ses-{ses}_{mprage}.nii.gz'
    sub_ses_segment = pjoin(bids, 'derivatives', 'SAMSEG', f'sub-{sub}', f'ses-{ses}', segment)
    sub_ses_transfo = pjoin(bids, 'derivatives', 'SAMSEG', f'sub-{sub}', f'ses-{ses}', transfo)
    
    # Perform Lesion Segmentation, Localisation and Volumetry        
    print('Start Preprocessing...')
    
    # Pre check
    ## Create corresponding directory if necessary
    os.makedirs(sub_ses_segment, exist_ok=True)
    os.makedirs(sub_ses_transfo, exist_ok=True)
    
    # Pre check
    ## check if preprocessing has already been done previously, and delete the results if necessary
    if os.path.isdir(pjoin(sub_ses_transfo, f'sub-{sub}_ses-{ses}_FLAIR_preprocessed')):
        os.remove(pjoin(sub_ses_transfo, f'sub-{sub}_ses-{ses}_FLAIR_preprocessed'))
    if pexists(pjoin(sub_ses_transfo, f'sub-{sub}_ses-{ses}_MPRAGE_preprocessed.nii.gz')):
        os.remove(pjoin(sub_ses_transfo, f'sub-{sub}_ses-{ses}_MPRAGE_preprocessed.nii.gz'))
          
    if not pexists(pjoin(sub_ses_anat, flair)):
        flair_bool = False
        print(f'No FLAIR for sub-{sub} ses-{ses}')
    if not pexists(pjoin(sub_ses_anat, mprage)):
        mprage_bool = False
        print(f'No MPRAGE for sub-{sub} ses-{ses}')
        
    if mprage_bool and flair_bool:
        shutil.copyfile(pjoin(sub_ses_anat, mprage), pjoin(sub_ses_transfo, f'sub-{sub}_ses-{ses}_MPRAGE.nii.gz'))
        
        try:
            print(f'Registration of FLAIR on MPRAGE for sub-{sub} ses-{ses}...')
            subprocess.Popen(f"docker run --rm --privileged -v {fs_license}:/usr/local/freesurfer/license.txt -v {sub_ses_anat}:/media/anat -v {sub_ses_transfo}:/media/transfo freesurfer:7.3.2 mri_coreg --mov {pjoin('/media/anat', flair)} --ref {pjoin('/media/anat', mprage)} --reg {pjoin('/media/transfo', 'flairToT1.lta')}", shell=True).wait()
            subprocess.Popen(f"docker run --rm --privileged -v {fs_license}:/usr/local/freesurfer/license.txt -v {sub_ses_anat}:/media/anat -v {sub_ses_transfo}:/media/transfo freesurfer:7.3.2 mri_vol2vol --mov {pjoin('/media/anat', flair)} --reg {pjoin('/media/transfo', 'flairToT1.lta')} --o {pjoin('/media/transfo', f'sub-{sub}_ses-{ses}_FLAIR.nii.gz')} --targ {pjoin('/media/anat', mprage)}", shell=True).wait()                                
            
        except Exception as e:
            print(f'[ERROR] - SAMSEG | {e} during registration of FLAIR on MPRAGE for sub-{sub}_ses{ses}!')
            return
        
    elif flair_bool:
        shutil.copyfile(pjoin(sub_ses_anat, flair), pjoin(sub_ses_transfo, f'sub-{sub}_ses-{ses}_FLAIR.nii.gz'))
    elif mprage_bool:
        shutil.copyfile(pjoin(sub_ses_anat, mprage), pjoin(sub_ses_transfo, f'sub-{sub}_ses-{ses}_MPRAGE.nii.gz'))
    else:
        return
            
        print('End Preprocessing!')
    
    ## Step2 : Segmentation using SAMSEG to compute lesion probability mask
        
    print('Start Segmentating lesions with SAMSEG...')
    
    # Actions
    ## Run SAMSEG with FLAIR and MPRAGE
    if mprage_bool and flair_bool:
        
        print('Segmenting with FLAIR & mprage')
        
        # Actions
        try:
            print(f'Running Samseg for sub-{sub} ses-{ses}...')
            
            subprocess.Popen(f"docker run --rm --privileged -v {fs_license}:/usr/local/freesurfer/license.txt -v {sub_ses_anat}:/media/anat -v {sub_ses_transfo}:/media/transfo freesurfer:7.3.2 mri_convert {pjoin('/media/transfo', f'sub-{sub}_ses-{ses}_FLAIR.nii.gz')} {pjoin('/media/transfo', f'sub-{sub}_ses-{ses}_FLAIR.mgz')}", shell=True).wait()
            
            subprocess.Popen(f"docker run --rm --privileged -v {fs_license}:/usr/local/freesurfer/license.txt -v {sub_ses_anat}:/media/anat -v {sub_ses_transfo}:/media/transfo freesurfer:7.3.2 mri_convert {pjoin('/media/transfo', f'sub-{sub}_ses-{ses}_MPRAGE.nii.gz')} {pjoin('/media/transfo', f'sub-{sub}_ses-{ses}_MPRAGE.mgz')}", shell=True).wait()
            
            subprocess.Popen(f"docker run --rm --privileged -v {fs_license}:/usr/local/freesurfer/license.txt -v {sub_ses_anat}:/media/anat -v {sub_ses_transfo}:/media/transfo -v {sub_ses_segment}:/media/segment freesurfer:7.3.2 run_samseg -i {pjoin('/media/transfo', f'sub-{sub}_ses-{ses}_MPRAGE.mgz')} -i {pjoin('/media/transfo', f'sub-{sub}_ses-{ses}_FLAIR.mgz')} --lesion --lesion-mask-pattern 0 1 --threads {cpus} -o {pjoin('/media/segment','SAMSEG_results')} --save-posteriors", shell=True).wait()
            
            subprocess.Popen(f"docker run --rm --privileged -v {fs_license}:/usr/local/freesurfer/license.txt -v {sub_ses_anat}:/media/anat -v {sub_ses_segment}:/media/segment freesurfer:7.3.2 mri_convert {pjoin('/media/segment', 'SAMSEG_results', 'posteriors', 'Lesions.mgz')} {pjoin('/media/segment', f'sub-{sub}_ses-{ses}_lesions.nii.gz')}", shell=True).wait()
            
        except Exception as e:
            print(f'[ERROR] - SAMSEG | {e} while running Samseg for sub-{sub}_ses{ses}!')
            return
        
    ## Run SAMSEG with only FLAIR 
    elif flair_bool:
        
        print('Segmenting with only FLAIR')
       
        # Actions
        try:
            print(f'Running Samseg for sub-{sub} ses-{ses}...')

            subprocess.Popen(f"docker run --rm --privileged -v {fs_license}:/usr/local/freesurfer/license.txt -v {sub_ses_anat}:/media/anat -v {sub_ses_transfo}:/media/transfo freesurfer:7.3.2 mri_convert {pjoin('/media/transfo', f'sub-{sub}_ses-{ses}_FLAIR.nii.gz')} {pjoin('/media/transfo', f'sub-{sub}_ses-{ses}_FLAIR.mgz')}", shell=True).wait()
            
            subprocess.Popen(f"docker run --rm --privileged -v {fs_license}:/usr/local/freesurfer/license.txt -v {sub_ses_anat}:/media/anat -v {sub_ses_transfo}:/media/transfo -v {sub_ses_segment}:/media/segment freesurfer:7.3.2 run_samseg -i {pjoin('/media/transfo', f'sub-{sub}_ses-{ses}_FLAIR.mgz')} --lesion --lesion-mask-pattern 0 1 --threads {cpus} -o {pjoin('/media/segment','SAMSEG_results')} --save-posteriors", shell=True).wait()
            
            subprocess.Popen(f"docker run --rm --privileged -v {fs_license}:/usr/local/freesurfer/license.txt -v {sub_ses_anat}:/media/anat -v {sub_ses_segment}:/media/segment freesurfer:7.3.2 mri_convert {pjoin('/media/segment', 'SAMSEG_results', 'posteriors', 'Lesions.mgz')} {pjoin('/media/segment', f'sub-{sub}_ses-{ses}_lesions.nii.gz')}", shell=True).wait()
            
        except Exception as e:   
            print(f'[ERROR] - SAMSEG | {e} while running Samseg for sub-{sub}_ses{ses}!')
            return
        
    elif mprage_bool:
        print('Segmenting with only MPRAGE')
        
        # Actions
        try:
            print(f'Running Samseg for sub-{sub} ses-{ses}...')

            subprocess.Popen(f"docker run --rm --privileged -v {fs_license}:/usr/local/freesurfer/license.txt -v {sub_ses_anat}:/media/anat -v {sub_ses_transfo}:/media/transfo freesurfer:7.3.2 mri_convert {pjoin('/media/transfo', f'sub-{sub}_ses-{ses}_MPRAGE.nii.gz')} {pjoin('/media/transfo', f'sub-{sub}_ses-{ses}_MPRAGE.mgz')}", shell=True).wait()
            
            subprocess.Popen(f"docker run --rm --privileged -v {fs_license}:/usr/local/freesurfer/license.txt -v {sub_ses_anat}:/media/anat -v {sub_ses_transfo}:/media/transfo -v {sub_ses_segment}:/media/segment freesurfer:7.3.2 run_samseg -i {pjoin('/media/transfo', f'sub-{sub}_ses-{ses}_MPRAGE.mgz')} --lesion --lesion-mask-pattern 0 1 --threads {cpus} -o {pjoin('/media/segment','SAMSEG_results')} --save-posteriors", shell=True).wait()
            
            subprocess.Popen(f"docker run --rm --privileged -v {fs_license}:/usr/local/freesurfer/license.txt -v {sub_ses_anat}:/media/anat -v {sub_ses_segment}:/media/segment freesurfer:7.3.2 mri_convert {pjoin('/media/segment', 'SAMSEG_results', 'posteriors', 'Lesions.mgz')} {pjoin('/media/segment', f'sub-{sub}_ses-{ses}_lesions.nii.gz')}", shell=True).wait()
            
        except Exception as e:   
            print(f'[ERROR] - SAMSEG | {e} while running Samseg for sub-{sub}_ses{ses}!')
            return
        
    ## Binarizing the lesion probability  mask    
    # Pre check
    ## check if files exist
    if not pexists(pjoin(sub_ses_segment, f'sub-{sub}_ses-{ses}_lesions.nii.gz')):
        file = pjoin(sub_ses_segment, f'sub-{sub}_ses-{ses}_lesions.nii.gz')
        print(f'[ERROR] - SAMSEG | FileNotFound: File "{file}" not Found !')
        return
    
    try:                
        # Actions
        print(f'Binarizing lesion probability mask for sub-{sub} ses-{ses}...')
       
        threshold = 0.5
        image = nib.load(pjoin(sub_ses_segment, f'sub-{sub}_ses-{ses}_lesions.nii.gz'))
        lesions = image.get_fdata()
        lesions[lesions >= threshold] = 1
        lesions[lesions < threshold] = 0

        
        lesions_nifti = nib.Nifti1Image(lesions, affine=image.affine)
        
        nib.save(lesions_nifti, pjoin(sub_ses_segment, f'sub-{sub}_ses-{ses}_lesions_binary.nii.gz'))

        
    except Exception as e:
        print(f'[ERROR] - SAMSEG | {e} while binarizing lesion probability mask for sub-{sub}_ses{ses}!')
        return
    
    print('End SAMSEG!')


def get_session_list(bids, subj, ses_details):
    """Helper function to get the list of sessions for a given subject."""
    sess = []
    if ses_details == 'all':
        for d in os.listdir(pjoin(bids, f'sub-{subj}')):
            if d.startswith('ses-'):
                sess.append(d.split('-')[1])
    else:
        for s in ses_details.split(','):
            if '-' in s:
                s0, s1 = map(int, s.split('-'))
                for si in range(s0, s1 + 1):
                    si_str = str(si).zfill(2)
                    if os.path.isdir(pjoin(bids, f'sub-{subj}', f'ses-{si_str}')):
                        sess.append(si_str)
            else:
                if os.path.isdir(pjoin(bids, f'sub-{subj}', f'ses-{s}')):
                    sess.append(s)
    return sess

def process_subject_range(bids, sub_range, ses_details):
    """Helper function to process a range of subjects."""
    subjects_and_sessions = []
    sub0, sub1 = map(int, sub_range.split('-'))
    for subi in range(sub0, sub1 + 1):
        subi_str = str(subi).zfill(3)
        if not os.path.isdir(pjoin(bids, f'sub-{subi_str}')):
            continue
        sess = get_session_list(bids, subi_str, ses_details)
        subjects_and_sessions.append((subi_str, sess))
    return subjects_and_sessions

def find_subjects_and_sessions(bids, sub, ses):
    subjects_and_sessions = []

    if sub == 'all':
        # Process all subjects
        for dirs in os.listdir(bids):
            if dirs.startswith('sub-'):
                subj = dirs.split('-')[1]
                sess = get_session_list(bids, subj, ses)
                subjects_and_sessions.append((subj, sess))
    else:
        # Process specified subjects
        for sub_item in sub.split(','):
            if '-' in sub_item:
                subjects_and_sessions.extend(process_subject_range(bids, sub_item, ses))
            else:
                if not os.path.isdir(pjoin(bids, f'sub-{sub_item}')):
                    continue
                sess = get_session_list(bids, sub_item, ses)
                subjects_and_sessions.append((sub_item, sess))
    
    return sorted(subjects_and_sessions)
    


if __name__ == '__main__':
    
    description = '''
bids_SAMSEG:
    This is a template python script to run SAMSEG
    '''
    
    usage = '\npython %(prog)s bids sub ses [OPTIONS]'
    
    parser = argparse.ArgumentParser(description=description, usage=usage)
    
    parser.add_argument('bids', type=str, help='path towards a bids formatted database')
    parser.add_argument('sub', type=str, help='sub ID or list of sub ID to process (e.g. 001,002). The keyword "all" will select all subjects of the database, while "-" allow to select subject ID in between two border (e.g. 001-010)')
    parser.add_argument('ses', type=str, help='ses ID or list of ses ID to process (e.g. 01,02). The keyword "all" will select all sessions of the database, while "-" allow to select session ID in between two border (e.g. 01-10)')
    parser.add_argument('--mprage', '-m', dest='mprage', type=str, help='name of the MPRAGE sequence', default=None, required=False)
    parser.add_argument('--flair', '-f', dest='flair', type=str, help='name of the FLAIR sequence', default=None, required=False)
    parser.add_argument('--cpus', dest='cpus', type=int, help='number of threads to run SAMSEG', default=8, required=False)
    parser.add_argument('--use-docker', dest='use_docker', help='use docker version of FreeSurfer instead of local one (default=False)', action='store_const', const=True, default=False, required=False)
    
    # Parse the arguments
    try:
        args = parser.parse_args()
    except SystemExit as e:
        # This block catches the SystemExit exception raised by argparse when required args are missing
        if e.code != 0:  # Non-zero code indicates an error
            parser.print_help()
        sys.exit(e.code)
        
    bids = args.bids
    
    subjects_and_sessions = find_subjects_and_sessions(bids, args.sub, args.ses)
    
    for sub, sess in subjects_and_sessions:
        for ses in sess:
            print(sub, ses)
            
            if args.use_docker:
                bids_SAMSEG_docker(bids, sub, ses, mprage=args.mprage, flair=args.flair, cpus=args.cpus)
            else:            
                bids_SAMSEG(bids, sub, ses, mprage=args.mprage, flair=args.flair, cpus=args.cpus)
    
    