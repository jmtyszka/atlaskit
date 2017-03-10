#!/usr/bin/env python3
"""
Create a report of intra and inter-observer atlas label statistics

Expects a label directory organized as follows:
<label_dir>/
  <observer A>/
    <template 1>
  <observer B>/


Usage
----
atlas.py -d <observer labels directory>
atlas.py -h

Authors
----
Mike Tyszka, Caltech Brain Imaging Center

Dates
----
2017-02-14 JMT From scratch

License
----
This file is part of atlaskit.

    atlaskit is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    atlaskit is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with atlaskit.  If not, see <http://www.gnu.org/licenses/>.

Copyright
----
2017 California Institute of Technology.
"""

__version__ = '0.2.0'

import os
import sys
import csv
import argparse
import nibabel as nib
import numpy as np
import pandas as pd
import multiprocessing as mp
import shutil
from glob import glob


def main():

    print()
    print('-------------------------------------------------------')
    print('Probablistic Atlas Construction with Similarity Metrics')
    print('-------------------------------------------------------')

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Calculate label similarity metrics for multiple observers')
    parser.add_argument('-d','--labeldir', help='Directory containing observer label subdirectories ["."]')
    parser.add_argument('-a','--atlasdir', help='Output atlas directory ["<labeldir>/atlas"]')
    parser.add_argument('-k','--key', help='ITK-SNAP label key text file ["<labeldir>/labels.txt"]')
    parser.add_argument('-l','--labels', required=False, type=parse_range, help='List of label indices to process (eg 1-5, 7-9, 12)')

    # Parse command line arguments
    args = parser.parse_args()

    if args.labeldir:
        label_dir = args.labeldir
        if not os.path.isdir(label_dir):
            print('Label directory does not exist (%s) - exiting' % label_dir)
            sys.exit(1)
    else:
        label_dir = os.path.realpath(os.getcwd())

    if args.atlasdir:
        atlas_dir = args.atlasdir
    else:
        atlas_dir = os.path.join(label_dir, 'atlas')

    # Check for label key file
    if os.path.isfile(args.key):
        label_keyfile = args.key
    else:
        print('* ITK-SNAP label key is missing (%s) - exiting' % args.key)
        sys.exit(1)

    # Safely create atlas directory
    if not os.path.isdir(atlas_dir):
        os.mkdir(atlas_dir)

    print('Label directory  : %s' % label_dir)
    print('Atlas directory  : %s' % atlas_dir)
    print('Label key file   : %s' % label_keyfile)

    # Save a copy of label key file in the atlas directory
    label_keyfile_save = os.path.join(atlas_dir, 'labels.txt')
    shutil.copyfile(label_keyfile, label_keyfile_save)

    # Load the label key as a data frame
    label_key = load_key(label_keyfile_save)

    # Init grand lists
    grand_labels = []
    vox_mm = []  # Voxel dimensions in mm
    vox_ul = []  # Voxel volume in mm^3 (microliters)
    affine_tx = []  # Nifti affine transform

    # Similarity metrics output files
    inter_metrics_csv = os.path.join(atlas_dir, 'inter_observer_metrics.csv')
    intra_metrics_csv = os.path.join(atlas_dir, 'intra_observer_metrics.csv')

    # Loop over observer directories
    # Any subdirectory of the label directory begining with "obs-"
    for obs_dir in glob(os.path.join(label_dir, "obs-*")):

        if os.path.isdir(obs_dir):

            print('Loading label images from %s' % obs_dir)

            # Init label image list for this observer
            obs_labels = []

            # Loop over all template label images
            for im in glob(os.path.join(obs_dir, '*.nii.gz')):

                # Load label image and add to list
                this_nii = nib.load(im)
                obs_labels.append(this_nii.get_data())

                # Save voxel dimensions, volume
                d = np.array(this_nii.header.get_zooms())
                vox_mm.append(d)
                vox_ul.append(d.prod())
                affine_tx.append(this_nii.get_affine())

            if len(obs_labels) > 0:
                print("  Loaded %d label images" % len(obs_labels))
                grand_labels.append(obs_labels)
            else:
                print("* No label images detected - skipping")

    # Voxel dimensions and volumes
    vox_mm, vox_ul = np.array(vox_mm), np.array(vox_ul)

    if not vox_mm.any():
        print("* No label images detected in %s - exiting" % label_dir)
        sys.exit(1)

    # Check for any variation in dimensions across templates and observers
    if any(np.nonzero(np.std(vox_mm, axis=1))):
        print('* Not all images have the same voxel dimensions - exiting')
        sys.exit(1)
    else:
        # Use dimensions from first image
        vox_mm = vox_mm[0]
        vox_ul = vox_ul[0]

    # Convert grand list to numpy array
    # -> labels[observer][template][x][y][z]
    print('Preparing labels')
    labels = np.array(grand_labels)

    # Limited list of labels to process
    if args.labels:
        label_nos = args.labels
    else:
        label_nos = np.int32(np.unique(labels))
        label_nos = np.delete(label_nos, np.where(label_nos == 0))  # Remove background label

    # Remove labels not present in key
    label_unknown = []
    for ll, label_no in enumerate(label_nos):
        if get_label_name(label_no, label_key) == 'Unknown':
            print('* Label %d unknown - removing from list' % label_no)
            label_unknown.append(ll)
    label_nos = np.delete(label_nos, label_unknown)

    # Count remaining labels
    n = len(label_nos)
    print('  Analyzing %d unique labels (excluding background)' % n)

    # Construct the probabilistic atlas from all labeled volumes
    prob_atlas = make_prob_atlas(labels, label_nos)

    # Save probabilistic atlas in label directory
    save_prob_atlas(prob_atlas, os.path.join(atlas_dir,'prob_atlas.nii.gz'), affine_tx[0])

    intra_metrics_all = []
    inter_metrics_all = []

    # Loop over each unique label value
    for label_no in label_nos:

        print('Analyzing label index %d' % label_no)

        # Current label mask
        label_mask = (labels == label_no)

        # Intra-observer metrics
        intra_metrics_all.append(intra_observer_metrics(label_mask, vox_mm))

        # Inter-observer metrics
        inter_metrics_all.append(inter_observer_metrics(label_mask, vox_mm))

    # Write metrics to report directory as CSV
    save_intra_metrics(intra_metrics_csv, intra_metrics_all, label_nos, label_key)
    save_inter_metrics(inter_metrics_csv, inter_metrics_all, label_nos, label_key)

    # Clean exit
    sys.exit(0)


def make_prob_atlas(labels, unique_labels):
    """

    Parameters
    ----------
    labels
    unique_labels

    Returns
    -------

    """

    print('Constructing probablistic atlas')

    # Get dimensions of label data
    n_obs, n_tmp, nx, ny, nz = labels.shape

    # Number of unique labels
    n = len(unique_labels)

    print('  Initializing probabilistic atlas')
    prob_atlas = np.zeros((nx, ny, nz, n), dtype='float32')

    # loop over each unique label value
    for m, label in enumerate(unique_labels):
        print('    Adding label %d' % label)
        label_mask = (labels == label).reshape(n_obs * n_tmp, nx, ny, nz)
        prob_atlas[:,:,:,m] = np.mean(label_mask, axis=0)

    return prob_atlas


def save_prob_atlas(prob_atlas, prob_atlas_fname, affine_tx):
    """

    Parameters
    ----------
    prob_atlas
    prob_atlas_fname
    affine_tx

    Returns
    -------

    """

    print('Saving probabilistic atlas to %s' % prob_atlas_fname)
    prob_nii = nib.Nifti1Image(prob_atlas, affine_tx)
    prob_nii.to_filename(prob_atlas_fname)


def intra_observer_metrics(label_mask, vox_mm):
    """
    Calculate within-observer Dice, Hausdorf and related metrics

    Parameters
    ----------
    label_mask: 5D numpy boolean array [observer][template][x][y][z]
    vox_mm: voxel dimensions in mm

    Returns
    -------
    intra_metrics: nobs x ntmp x ntmp nested list
    """

    # Dimensions
    nobs, ntmp, nx, ny, nz = label_mask.shape

    # Init grand intra-observer metrics list
    intra_metrics = []

    print('  Calculating intra-observer similarity metrics :', end='')

    for obs in range(0,nobs):

        print(' %d' % obs, end='')

        # Results list for current observer
        obs_res = []

        for ta in range(0, ntmp):

            mask_a = label_mask[obs, ta, :, :, :]
            data_a = []

            for tb in range(0, ntmp):

                mask_b = label_mask[obs,tb,:,:,:]
                data_a.append((mask_a, mask_b, vox_mm))

            # Run similarity metric function in parallel on template A data list
            with mp.Pool(mp.cpu_count()-2) as pool:
                res = pool.starmap(similarity, data_a)

            # Add to current observer results
            obs_res.append(res)

        # Add observer results to grand list
        intra_metrics.append(obs_res)

    print()

    return intra_metrics


def inter_observer_metrics(label_mask, vox_mm):
    """
     Calculate between-observer Dice, Hausdorf and related metrics

     Parameters
     ----------
     label_mask: 5D numpy boolean array [observer][template][x][y][z]
     vox_mm: voxel dimensions in mm

     Returns
     -------
     inter_metrics: ntmp x nobs x nobs nested list
     """

    # Dimensions
    nobs, ntmp, nx, ny, nz = label_mask.shape

    # Init grand inter-observer metrics list
    inter_metrics = []

    print('  Calculating inter-observer similarity metrics :', end='')

    for tmp in range(0, ntmp):

        print(' %d' % tmp, end='')

        # Results list for this template
        tmp_res = []

        for obs_a in range(0, nobs):

            mask_a = label_mask[obs_a, tmp, :, :, :]

            data_a = []

            for obs_b in range(0, nobs):

                mask_b = label_mask[obs_b, tmp, :, :, :]

                data_a.append((mask_a, mask_b, vox_mm))

            # Run similarity metric function in parallel on data list
            with mp.Pool(mp.cpu_count()-2) as pool:
                res = pool.starmap(similarity, data_a)

            # Add to current template results
            tmp_res.append(res)

        # Add template results to grand list
        inter_metrics.append(tmp_res)

    print()

    return inter_metrics


def save_intra_metrics(fname, intra_metrics, label_nos, label_key):
    """

    Parameters
    ----------
    fname: CSV filename
    intra_metrics: nobs x ntmp x ntmp nested list
    label_nos:
    label_key:

    Returns
    -------

    """

    print('Saving intra-observer metrics to %s' % fname)

    # Preferred method for safe opening CSV file in Python 3
    with open(fname, "w", newline='') as f:

        writer = csv.writer(f)

        # Column headers
        writer.writerow(('labelName','labelNo','observer','tmpA','tmpB','dice','hausdorff','nA','nB'))

        for idx, m_idx in enumerate(intra_metrics):
            label_no = label_nos[idx]
            label_name = get_label_name(label_no, label_key)
            for obs, m_obs in enumerate(m_idx):
                for tA, m_ta in enumerate(m_obs):
                    for tB, m_tb in enumerate(m_ta):
                        writer.writerow((label_name, label_no, obs, tA, tB) + m_tb)


def save_inter_metrics(fname, inter_metrics, label_nos, label_key):
    """

    Parameters
    ----------
    fname: CSV filename
    inter_metrics: ntmp x nobs x nobs nested list
    label_nos:
    label_key:

    Returns
    -------

    """

    print('Saving inter-observer metrics to %s' % fname)

    with open(fname, "w", newline='') as f:

        writer = csv.writer(f)

        # Column headers
        writer.writerow(('labelName','labelNo','template','obsA','obsB','dice','hausdorff','nA','nB'))

        for idx, m_idx in enumerate(inter_metrics):
            label_no = label_nos[idx]
            label_name = get_label_name(label_no, label_key)
            for tmp, m_tmp in enumerate(m_idx):
                for obsA, m_oa in enumerate(m_tmp):
                    for obsB, m_ob in enumerate(m_oa):
                        writer.writerow((label_name, label_no, tmp, obsA, obsB) + m_ob)


def similarity(mask_a, mask_b, vox_mm):

    # Count voxels in each mask
    na, nb = np.sum(mask_a), np.sum(mask_b)

    # Only calculate stats if labels present in A or B
    if na > 0 or nb > 0:

        # Find intersection and union of A and B masks
        a_and_b = np.logical_and(mask_a, mask_b)
        a_or_b = np.logical_or(mask_a, mask_b)

        # Count voxels in intersection and union
        n_a_and_b, n_a_or_b = np.sum(a_and_b), np.sum(a_or_b)

        # Similarity metrics
        dice = 2.0 * n_a_and_b / float(na + nb)
        haus = hausdorff_distance(mask_a, mask_b, vox_mm)
    else:
        dice, haus = np.nan, np.nan

    return dice, haus, na, nb


def hausdorff_distance(A, B, vox_mm):
    """
    Calculate the hausdorff_distance distance in mm between two binary masks in 3D

    Parameters
    ----------
    A : 3D numpy array
        Binary mask A
    B : 3D numpy array
        Binary mask B
    vox_mm : numpy array
        voxel dimensions in mm

    Returns
    -------
    H : float
        hausdorff_distance distance between labels
    """

    # Create lists of all True points in both masks
    xA, yA, zA = np.nonzero(A)
    xB, yB, zB = np.nonzero(B)

    # Count elements in each point set
    nA = xA.size
    nB = xB.size

    if nA > 0 and nB > 0:

        # Init min dr to -1 for all points in A
        min_dr = -1.0 * np.ones([nA])

        for ac in range(0,nA):

            dx = (xA[ac] - xB[:]) * vox_mm[0]
            dy = (yA[ac] - yB[:]) * vox_mm[1]
            dz = (zA[ac] - zB[:]) * vox_mm[2]
            min_dr[ac] = np.min(np.sqrt(dx**2 + dy**2 + dz**2))

        # Find maximum over A of the minimum distances A to B
        H = np.max(min_dr)

    else:

        H = np.nan

    return H


def get_template_ids(label_dir, obs):

    obs_dir = os.path.join(label_dir, obs)

    if os.path.isdir(obs_dir):

        ims = glob(os.path_join(obs_dir, '*.nii.gz'))

        print(ims)


def parse_range(astr):
    '''
    Parse compound list of integers and integer ranges

    Parameters
    ----------
    astr

    Returns
    -------

    '''
    result = set()
    for part in astr.split(','):
        x = part.split('-')
        result.update(range(int(x[0]), int(x[-1]) + 1))

    return sorted(result)


def load_key(key_fname):
    """
    Parse an ITK-SNAP label key file

    Parameters
    ----------
    key_fname: ITK-SNAP label key filename

    Returns
    -------
    key: Data table containing ITK-SNAP style label key
    """

    # Import key as a data table
    # Note the partially undocumented delim_whitespace flag
    key = pd.read_table(key_fname,
                         comment='#',
                         header=None,
                         names=['Index','R','G','B','A','Vis','Mesh','Name'],
                         delim_whitespace=True)

    return key


def get_label_name(label_idx, label_key):

    label_name = 'Unknown'

    for i, idx in enumerate(label_key.Index):
        if label_idx == idx:
            label_name = label_key.Name[i]

    return label_name


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
