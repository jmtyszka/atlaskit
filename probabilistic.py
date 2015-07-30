#!/opt/local/bin/python
"""
Construct a probabilistic atlas from a set of N label images each containing
M unique labels. The final atlas will be a 4D float image (nx x ny x nz x M)
scaled from 0.0 to 1.0

Usage
----
probabilistic.py -o <Output filename> <List of N label image volumes>
probabilistic.py -h

Example
----
>>> probabilistic.py Labels_1.nii.gz Labels_2.nii.gz ... Labels_N.nii.gz

Authors
----
Mike Tyszka, Caltech Brain Imaging Center

Dates
----
2015-07-29 JMT From scratch

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


def main():
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Construct probabilistic atlas from label volumes')
    parser.add_argument('-o', '--output', help='Output atlas filename')
    parser.add_argument('label_files', nargs='+', help='Space-separated list of label filenames')

    # Parse command line arguments
    args = parser.parse_args()
    
    if args.output:
        prob_file = args.output
    else:
        prob_file = 'Probabilistic_Atlas.nii.gz'
        
    label_files = args.label_files
    
    # Count number of label files
    N = len(label_files)
    
    for i, fname in enumerate(label_files):
        
        print('  Adding label volume ' + fname)
        
        # Open current label volume
        label_nii = nib.load(fname)

        # Get data from current label volume
        labels = label_nii.get_data()

        # Init from first volume       
        if i < 1:
            
            print('  Initializing probabilistic atlas')

            # Construct list of unique label values in first volume
            # NOTE: should be identical list for all volumes
            unique_labels = np.unique(labels)
            
            # Remove background label (0) from list
            unique_labels = np.delete(unique_labels, np.where(unique_labels == 0))
            
            # Count number of labels
            M = len(unique_labels)
            
            print('  Identified %d unique labels' % M)
            
            # Make space for final 4D prob atlas
            nx,ny,nz = label_nii.header.get_data_shape()
            prob = np.zeros((nx,ny,nz,M), dtype='float32')
            
            # Grab affine transform from first label volume
            T = label_nii.get_affine()
    
        # loop over each unique label value
        for m, label in enumerate(unique_labels):
            
            # Create mask for current label
            label_mask = (labels == label)
                
            # Add mask to 4D
            prob[:,:,:,m] += label_mask
    
    # Normalize probabilities to [0,1]
    prob /= float(N)
    
    # Write 4D probabilistic atlas
    print('Saving probabilistic atlas to %s' % prob_file)
    prob_nii = nib.Nifti1Image(prob, T)
    prob_nii.to_filename(prob_file)
    
    # Clean exit
    sys.exit(0)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
