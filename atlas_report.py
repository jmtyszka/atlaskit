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
atlas_report.py -d <observer labels directory>
atlas_report.py -h

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
from glob import glob


def main():

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Report intra and inter-observer atlas label statistics')
    parser.add_argument('-d','--labeldir', required=True, help='Directory containing observer labels')
    parser.add_argument('-k', '--key', required=False, help='ITK-SNAP label key file [Atlas_Labels.txt]')

    # Parse command line arguments
    args = parser.parse_args()

    label_dir = args.labeldir

    # Check for label directory existence
    if not os.path.isdir(label_dir):
        print('Label directory does not exist (%s) - existing' % label_dir)
        sys.exit(1)

    # Load and parse label key if provided
    if args.key:
        label_key = load_key(os.path.join(label_dir, args.key))
    else:
        label_key = []

    # Get list of observers
    observers = next(os.walk(label_dir))[1]

    # Init grand lists
    grand_labels = []
    vox_mm = []  # Voxel dimensions in mm
    vox_ul = []  # Voxel volume in mm^3 (microliters)

    # Create report directory alongside label directory
    report_dir = label_dir + '.report'
    if not os.path.isdir(report_dir):
        os.mkdir(report_dir)
    print('Report directory : %s' % report_dir)

    # Inter/intra-observer metrics CSV files
    inter_metrics_csv = os.path.join(report_dir, 'inter_observer_metrics.csv')
    intra_metrics_csv = os.path.join(report_dir, 'intra_observer_metrics.csv')

    # Compute metrics if required, otherwise skip to next section
    if not os.path.isfile(inter_metrics_csv) and not os.path.isfile(intra_metrics_csv):

        # Loop through all subdirectories of label directory
        for obs in observers:

            print('Loading label images for observer %s ... ' % obs, end='')

            # Init label image list for this observer
            obs_labels = []

            # Observer subdirectory
            obs_dir = os.path.join(label_dir, obs)

            # Loop over all template label images
            for im in glob(os.path.join(obs_dir, '*.nii.gz')):

                # Load label image and add to list
                this_nii = nib.load(im)
                obs_labels.append(this_nii.get_data())

                # Save voxel dimensions, volume
                d = np.array(this_nii.header.get_zooms())
                vox_mm.append(d)
                vox_ul.append(d.prod())

            print("%d label images loaded" % len(obs_labels))

            # Add observer labels to grand labels
            grand_labels.append(obs_labels)

        # Voxel dimensions and volumes
        vox_mm, vox_ul = np.array(vox_mm), np.array(vox_ul)

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

        # Find unique label numbers within data
        print('Determining unique label indices')
        unique_labels = np.unique(labels)

        print('Found %d unique label indices' % len(unique_labels))

        intra_metrics_all = []
        inter_metrics_all = []

        # loop over each unique label value
        # for label_idx in unique_labels:
        for label_idx in range(0,3):

            if label_idx > 0:

                # Find label name if provided
                if np.any(label_key):
                    label_name = get_label_name(label_idx, label_key)
                else:
                    label_name = 'Unknown'

                print('Analysing label %d (%s)' % (label_idx, label_name))

                # Current label mask
                label_mask = (labels == label_idx)

                # Intra-observer metrics
                intra_metrics_all.append((intra_observer_metrics(label_mask, vox_mm), label_idx, label_name))

                # Inter-observer metrics
                inter_metrics_all.append((inter_observer_metrics(label_mask, vox_mm), label_idx, label_name))

        # Write metrics to report directory as CSV
        save_intra_metrics(intra_metrics_csv, intra_metrics_all)
        save_inter_metrics(inter_metrics_csv, inter_metrics_all)


    # Clean exit
    sys.exit(0)


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


def save_intra_metrics(fname, intra_metrics):
    """

    Parameters
    ----------
    fname: CSV filename
    intra_metrics: nobs x ntmp x ntmp nested list

    Returns
    -------

    """

    with open(fname, "wb") as f:
        writer = csv.writer(f)

        for m_obs in intra_metrics:


def save_inter_metrics(fname, inter_metrics):
    """

    Parameters
    ----------
    fname: CSV filename
    inter_metrics: ntmp x nobs x nobs nested list

    Returns
    -------

    """

    with open(fname, "wb") as f:
        writer = csv.writer(f)


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


def load_key(key_fname):
    '''
    Parse an ITK-SNAP label key file

    Parameters
    ----------
    key_fname

    Returns
    -------

    '''

    # Import key as a data table
    # Note the partially undocumented delim_whitespace flag
    data = pd.read_table(key_fname,
                         comment='#',
                         header=None,
                         names=['Index','R','G','B','A','Vis','Mesh','Name'],
                         delim_whitespace=True)

    return data


def get_label_name(label_idx, label_key):
    '''
    Search label key for label index and return name

    Parameters
    ----------
    label_idx
    label_key

    Returns
    -------

    '''

    label_name = 'Unknown Label'

    for i, idx in enumerate(label_key.Index):
        if label_idx == idx:
            label_name = label_key.Name[i]

    return label_name


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


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
