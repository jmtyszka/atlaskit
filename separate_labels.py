#!/usr/bin/env python3
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
2015-07-29 JMT Speed up mask generation, use zero-padded output indexing

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

__version__ = '0.2.0'

import os
import sys
import argparse
import nibabel as nib
import numpy as np


def main():
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Smooth one or more atlas labels')
    parser.add_argument('in_file', help="source atlas labels filename")
    
    args = parser.parse_args()
    
    in_file = args.in_file
    
    # Convert relative to absolute path
    in_file = os.path.abspath(in_file)
    
    # Load the source atlas image
    print('Opening %s' % in_file)
    in_nii = nib.load(in_file)
    
    # Load label image
    src_labels = in_nii.get_data()

    # Construct list of unique label values in image
    unique_labels = np.unique(src_labels)

    # loop over each unique label value
    for label in unique_labels:
        
        if label > 0:
        
            # Create mask for current label value
            out_mask = (src_labels == label).astype(int)
            
            # Construct output filename. Use zero-padded indexing
            out_file = in_file.strip('.nii.gz') + '_' + '{0:04d}'.format(label) + '.nii.gz'
            
            # Save smoothed labels image
            print('Saving label %d to %s' % (label, out_file))
            out_nii = nib.Nifti1Image(out_mask, in_nii.get_affine())
            out_nii.to_filename(out_file)
    
    print('Done')
    
    # Clean exit
    sys.exit(0)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
