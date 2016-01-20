#!/opt/local/bin/python
"""
Calculate the local intensity gradient using a 3D Sobel filter

Usage
----
sobel.py -i <input image> -o <output image>
sobel.py -h

Example
----
>>> sobel.py -i T1.nii.gz -o T1_sobel.nii.gz

Authors
----
Mike Tyszka, Caltech Brain Imaging Center

Dates
----
2016-01-06 JMT From scratch 

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
2016 California Institute of Technology.
"""

__version__ = '0.1.0'

import sys
import argparse
import numpy as np
import nibabel as nib
from scipy.ndimage.filters import sobel


def main():
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Sobel filter 3D image')
    parser.add_argument('-i','--in_file', help="source image filename")
    parser.add_argument('-o','--out_file', help="Sobel filtered image filename")

    args = parser.parse_args()

    in_file = args.in_file
    out_file = args.out_file
        
    # Load the source image
    print('Opening %s' % in_file)
    in_nii = nib.load(in_file)
    
    # Load label image
    print('Loading image')
    src_img = in_nii.get_data()
    src_img = src_img.astype(float)
    
    # Smooth target label region
    print('  Sobel gradients')
    Sx = sobel(src_img, axis=0)
    Sy = sobel(src_img, axis=1)
    Sz = sobel(src_img, axis=2)
    
    # Take gradient magnitude
    print('  Sobel gradient magnitude')
    # out_img = np.sqrt(Sx**2 + Sy**2 + Sz**2)
    out_img = np.sqrt(Sx**2 + Sz**2)
    
    # Save Sobel image
    print('Saving Sobel image %s' % out_file)
    out_nii = nib.Nifti1Image(out_img, in_nii.get_affine())
    out_nii.to_filename(out_file)
    
    print('Done')
    
    # Clean exit
    sys.exit(0)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
