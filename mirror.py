#!/opt/local/bin/python
"""
Flip 3D (xyz) or 4D (xyzt) data volume in x, WITHOUT ALTERING THE HEADER
The usual caveats apply. ONLY RUN THIS ON ABSTRACTED DATA (eg atlases)

Usage
----
mirror.py <3D or 4D image>
mirror.py -h

Example
----
>>> mirror.py atlas_labels.py

Authors
----
Mike Tyszka, Caltech Brain Imaging Center

Dates
----
2015-07-31 JMT From scratch
2016-03-25 JMT Explicit output filename

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


def main():
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Mirror data in x without header adjustment')
    parser.add_argument('-i','--in_file', required=True, help="Input Nifti 3D or 4D volume")
    parser.add_argument('-o','--out_file', required=True, help="Output mirrored version of input image")

    # Parse command line arguments
    args = parser.parse_args()
    in_file = args.in_file
    out_file = args.out_file
    
    # Convert relative to absolute path
    in_file = os.path.abspath(in_file)
    out_file = os.path.abspath(out_file)

    # Open and load Nifti image
    in_nii = nib.load(in_file)
    in_data = in_nii.get_data()

    nd = in_data.ndim

    # Flip image data in first dimension (x)
    if nd == 3:
        out_data = in_data[::-1,:,:]
    elif nd == 4:
        out_data = in_data[::-1,:,:,:]

    # Write x-mirrored data with identical header
    print('Saving x-mirrored image to %s' % out_file)
    out_nii = nib.Nifti1Image(out_data, in_nii.get_affine())
    out_nii.to_filename(out_file)
    
    # Clean exit
    sys.exit(0)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
