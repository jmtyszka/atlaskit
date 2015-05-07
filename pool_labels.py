#!/opt/local/Library/Frameworks/Python.framework/Versions/2.7/Resources/Python.app/Contents/MacOS/Python
"""
Pool/Combine one or more labels within an integer atlas volume

Usage
----
pool_labels.py <input label image> <output label image> <input label numbers> <output label number>
pool_labels.py -h

Example
----
>>> pool_labels.py atlas.nii.gz atlas_1.nii.gz 1 12 13 14

Authors
----
Mike Tyszka, Caltech Brain Imaging Center
Wolfgang M. Pauli, Caltech
Kabir Brar, Caltech 

Dates
----
2015-05-05 WMP and KB from scratch

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
    parser = argparse.ArgumentParser(description='Pool one or more atlas labels')
    parser.add_argument('in_file', help="source atlas labels filename")
    parser.add_argument('out_file', help="pooled atlas labels filename")
    parser.add_argument('out_label', metavar='out_label', type=int, help='one label number to renumber in_labels to')
    parser.add_argument('in_labels', metavar='in_label', type=int, nargs='+',
                        help='two label numbers to change')

    args = parser.parse_args()

    in_file = args.in_file
    out_file = args.out_file
    in_labels = args.in_labels
    out_label = args.out_label
        
    # Load the source atlas image
    print('Opening %s' % in_file)
    in_nii = nib.load(in_file)
    
    # Load label image
    print('Loading labels')
    src_labels = in_nii.get_data()
        
    print('Pooling desired labels')
    out_labels = np.zeros_like(src_labels)
    for label in in_labels:
        # Extract target label as a boolean mask
        out_labels += (src_labels == label)

    # Renumber target label region
    out_labels *= out_label
  
    print('Copying all other labels from original')
    for label in np.unique(src_labels):
        if (label in in_labels) or label == out_label or label == 0:
            continue
        else:
            ind = np.where(src_labels == label)
            out_labels[ind] = label
    
    # Save changed labels image
    print('Saving changed labels to %s' % out_file)
    out_nii = nib.Nifti1Image(out_labels, in_nii.get_affine())
    out_nii.to_filename(out_file)
    
    print('Done')
    
    # Clean exit
    sys.exit(0)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
