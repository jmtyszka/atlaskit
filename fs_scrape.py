#!/usr/bin/env python3
"""
Scrape all cortical (aparc 2009) and subcortical (aseg) volumes into a single CSV

2019-05-08 Mike Tyszka 
"""

import os
import sys
import pandas as pd
from glob import glob


def load_fs_cort_vols(subj_dir, hemi='lh'):

    fname = os.path.join(subj_dir, 'stats', '{}.aparc.a2009s.stats'.format(hemi)) 

    if os.path.isfile(fname):

        print('Loading {}'.format(fname))

        vols = pd.read_csv(fname,
                           sep=" ",
                           header=None,
                           usecols=[0, 3],
                           skipinitialspace=True,
                           comment='#')        
        vols.columns = ['label', '{}_gm_vol_mm3'.format(hemi)]

    else:

        print('{} not found - skipping'.format(fname))
        vols = []

    return vols


def load_fs_subcort_vols(subj_dir):

    fname = os.path.join(subj_dir, 'stats', 'aseg.stats')

    if os.path.isfile(fname):

        print('Loading {}'.format(fname))

        vols = pd.read_csv(fname,
                           sep=" ",
                           header=None,
                           usecols=[3, 4],
                           skipinitialspace=True,
                           comment='#')
        vols.columns = ['vol_mm3', 'label']

    else:

        print('{} not found - skipping'.format(fname))
        vols = []

    return vols


# Key directories
bids_dir = os.path.abspath('..')
der_dir = os.path.join(bids_dir, 'derivatives')
fs_dir = os.path.join(der_dir, 'freesurfer')
vols_dir = os.path.join(der_dir, 'volumes')

if not os.path.isdir(vols_dir):
    os.makedirs(vols_dir, exist_ok=True)

vol_list = []

# Freesurfer subject results in derivatives/freesurfer/sub-<Subject ID>_core<Core Version>/
# <Core Version> = 1 or 2

for subj_dir in glob(os.path.join(fs_dir, 'sub-CC*')):

    subj_str = os.path.basename(subj_dir)
    subj_id, core_ver = subj_str.split('_')
    _, subj_id = subj_id.split('-')
    core_ver = core_ver[-1]

    print('Subject {} Core {}'.format(subj_id, core_ver))

    # Freesurfer cortical volumes file
    lh_cort_vols = load_fs_cort_vols(subj_dir, 'lh')
    rh_cort_vols = load_fs_cort_vols(subj_dir, 'rh')
    subcort_vols = load_fs_subcort_vols(subj_dir)


    

    

