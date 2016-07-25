#!/usr/bin/env python3
"""
Smooth one or more labels within an integer atlas volume

Usage
----
smooth_labels.py -i <input label image> -o <output label image> [label numbers]
smooth_labels.py -h

Example
----
>>> smooth_labels.py -i atlas.nii.gz -o atlas_smooth_5.nii.gz 5 10 11

Authors
----
Mike Tyszka, Caltech Brain Imaging Center

Dates
----
2015-04-07 JMT From scratch
2015-12-08 JMT Update command line arguments and port to python 3

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


def main():
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Smooth one or more atlas labels')
    parser.add_argument('-i','--in_file', help="source atlas labels filename")
    parser.add_argument('-o','--out_file', help="smoothed atlas labels filename")
    parser.add_argument('labels', metavar='N', type=int, nargs='+',
                        help='label numbers to smooth')

    args = parser.parse_args()

    in_file = args.in_file
    out_file = args.out_file
    labels = args.labels
        
    # Load the source atlas image
    print('Opening %s' % in_file)
    in_nii = nib.load(in_file)
    
    # Load label image
    print('Loading labels')
    src_labels = in_nii.get_data()
    
    # Duplicate into output image
    print('Creating new label image')
    out_labels = src_labels.copy()
    
    for label in labels:
        
        print('  Smoothing label %d' % label)
    
        # Extract target label as a boolean mask
        print('  Identifying target label region %d' % label)
        label_mask = (src_labels == label)
    
        # Smooth target label region
        print('  Gaussian smoothing original target label')
        label_mask_smooth = gaussian_filter(label_mask.astype(float), sigma=1.0)
    
        # Normalize smoothed intensities
        label_mask_smooth = label_mask_smooth / label_mask_smooth.max()
    
        # Threshold smoothed mask at 0.5 to create new boolean mask
        print('  Thresholding smoothed label')
        label_mask_smooth = label_mask_smooth > 0.5
    
        # Replace unsmoothed with smoothed label, overwriting other labels
        print('  Inserting smoothed target label')
        out_labels[label_mask] = 0
        out_labels[label_mask_smooth] = label
    
    # Save smoothed labels image
    print('Saving smoothed labels to %s' % out_file)
    out_nii = nib.Nifti1Image(out_labels, in_nii.get_affine())
    out_nii.to_filename(out_file)
    
    print('Done')
    
    # Clean exit
    sys.exit(0)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
