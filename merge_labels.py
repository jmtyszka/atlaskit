#!/usr/bin/python
"""
Merge separate segmentation volumens, with one ONLY one label per volume, into one.

Usage
----
merge_labels.py <out label image> <first input label image>  <second input label image> ...
merge_labels.py -h

Example
----
>>> merge_labels.py atlas.nii.gz seg_1.nii.gz seg_2.nii.gz 

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
    parser.add_argument('out_file', help="smoothed atlas labels filename")
    parser.add_argument('in_files', metavar='N', type=str, nargs='+',
                        help='label numbers to smooth')
    
    args = parser.parse_args()
    
    out_file = args.out_file
    in_files = args.in_files
    
    # Load label image
    print('Opening first input file to use as reference')
    in_nii = nib.load(in_files[0])
    out_labels = np.zeros_like(in_nii.get_data())

    for i in xrange(len(in_files)):
        in_file = in_files[i]

        # Load the source atlas image
        print('Processing %s' % in_file)
        in_nii = nib.load(in_file)
            
        # Load label image
        src_labels = in_nii.get_data()
        out_labels[np.where(src_labels)] = i + 1


    # Save smoothed labels image
    print('Saving merged labels to %s' % out_file)
    out_nii = nib.Nifti1Image(out_labels, in_nii.get_affine())
    out_nii.to_filename(out_file)

    
    print('Done')
    
    # Clean exit
    sys.exit(0)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
