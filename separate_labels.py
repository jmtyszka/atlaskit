#!/usr/bin/python

"""
Write separate volume for each label found in the input image

Usage
----
separate_labels.py <input label image>
separate_labels.py -h

Example
----
>>> separate_labels.py atlas.nii.gz 

Authors
----
Wolfgang M. Pauli, Caltech, Division of Humaninities and Social Sciences

Dates
----
2015-05-02 WMP From scratch

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
from scipy.ndimage.filters import gaussian_filter
import nibabel as nib
import numpy as np


def main():
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Smooth one or more atlas labels')
    parser.add_argument('in_file', help="source atlas labels filename")
    
    args = parser.parse_args()
    
    in_file = args.in_file
    
    # Load the source atlas image
    print('Opening %s' % in_file)
    in_nii = nib.load(in_file)
    
    # Load label image
    src_labels = in_nii.get_data()
    
    unique_labels = np.unique(src_labels)

    for label in unique_labels[1:len(unique_labels)]:
        out_labels = np.zeros_like(src_labels)

        ind = np.where(src_labels == label)
        
        out_labels[ind] = 1
        
        out_file = in_file.strip('.nii.gz') + '_' + str(int(label)) + '.nii.gz'
        
        # Save smoothed labels image
        print('Saving label %s to %s' % (str(int(label)), out_file))
        out_nii = nib.Nifti1Image(out_labels, in_nii.get_affine())
        out_nii.to_filename(out_file)
    
    print('Done')
    
    # Clean exit
    sys.exit(0)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
