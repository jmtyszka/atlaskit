#!/usr/bin/env python3
"""
Calculate generalized Hausdorff distance between corresponding labels in A and B

Usage
----
hausdorff.py <labelsA> <labelsB>
haussdorf.py -h

Example
----
>>> hausdorff.py labelsA.nii.gz labelsB.nii.gz

Authors
----
Mike Tyszka, Caltech Brain Imaging Center

Dates
----
2016-03-15 JMT Adapt from dice.py

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
2015 California Institute of Technology.
"""

__version__ = '0.1.0'

import sys
import argparse
import nibabel as nib
import numpy as np
import pandas as pd


def main():

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Hausdorff distance between corresponding labels in two volumes')
    parser.add_argument('-a','--labelsA', required=True, help='Labeled volume A')
    parser.add_argument('-b','--labelsB', required=True, help='Labeled volume B')
    parser.add_argument('-k','--labelsKey', help='ITK-SNAP label key [optional]')


    # Parse command line arguments
    args = parser.parse_args()

    labelsA = args.labelsA
    labelsB = args.labelsB

    if args.labelsKey:
        labelsKey = args.labelsKey
    else:
        labelsKey = ''

    # Load labeled volumes
    A_nii, B_nii = nib.load(labelsA), nib.load(labelsB)
    A_labels, B_labels = A_nii.get_data(), B_nii.get_data()

    # Load and parse label key
    label_key = LoadKey(labelsKey)

    # Voxel volume in mm^3 (microliters) (from volume A)
    atlas_vox_vol_ul = np.array(A_nii.header.get_zooms()).prod()

    # Colume headers
    print('%24s,%8s,%8s,%8s,%10s,%10s,%10s,%10s,' %
        ('Label', 'Index', 'nA', 'nB', 'vA_ul', 'vB_ul', 'Hausdorff'))

    # Construct list of unique label values in image
    unique_labels = np.unique(A_labels)

    # loop over each unique label value
    for label_idx in unique_labels:

        if label_idx > 0:

            # Find label name
            label_name = LabelName(label_idx, label_key)

            # Create label mask from A and B volumes
            A_mask = (A_labels == label_idx)
            B_mask = (B_labels == label_idx)

            # Count voxels in each mask
            nA, nB = np.sum(A_mask), np.sum(B_mask)

            # Only calculate stats if labels present in A or B
            if nA > 0 or nB > 0:

                # Similarity coefficients
                H = Hausdorff(A_mask, B_mask)

                # Absolute volumes of label in A and B
                A_vol_ul = np.sum(A_mask) * atlas_vox_vol_ul
                B_vol_ul = np.sum(B_mask) * atlas_vox_vol_ul

                print('%24s,%8d,%8d,%8d,%10.3f,%10.3f,%10.3f,%10.3f,' %
                    (label_name, label_idx, nA, nB, A_vol_ul, B_vol_ul, H))

    # Clean exit
    sys.exit(0)

def Hausdorff(A, B):
    """
    Calculate the Hausdorff distance between two binary masks in 3D

    Parameters
    ----------
    A : 3D numpy array
    B : 3D numpy array

    Returns
    -------
    H : float
        Hausdorff distance between labels
    """

    # Create lists of all True points in both masks


    return H

def LoadKey(key_fname):
    '''
    Parse an ITK-SNAP label key file
    '''

    # Import key as a data table
    # Note the partially undocumented delim_whitespace flag
    data = pd.read_table(key_fname,
                         comment='#',
                         header=None,
                         names=['Index','R','G','B','A','Vis','Mesh','Name'],
                         delim_whitespace=True)

    return data

def LabelName(label_idx, label_key):
    '''
    Search label key for label index and return name
    '''

    label_name = 'Unknown Label'

    for i, idx in enumerate(label_key.Index):
        if label_idx == idx:
            label_name = label_key.Name[i]

    return label_name


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
