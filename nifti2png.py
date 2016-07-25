#!/usr/bin/env python3
"""
Convert single Nifti-1 3D or 4D volume to set of PNG stacks.

Usage
----
nifti2png.py -i <Nifti filename> -o <PNG image stub>
nifti2png.py -h

Example
----
>>> nifti2png.py -i mynifti.nii.gz -o mypngs

Authors
----
Mike Tyszka, Caltech Brain Imaging Center

Dates
----
2016-06-28 JMT Adapt from nifti2jpg.py

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

import os
import sys
import cv2
import argparse
import nibabel as nib
import numpy as np


def main():

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Convert 3D Nifti volume to JPEG stack')
    parser.add_argument('-i', '--nii_file', required=True, help='3D or 4D Nifti image')
    parser.add_argument('-o', '--png_stub', required=False, help='PNG stack stub')

    # Parse arguments
    args = parser.parse_args()
    nii_file = args.nii_file

    if args.png_stub:
        png_stub = args.png_stub
    else:
        png_stub = "Images"

    # Strip extension from Nifti filename for use as a stub
    nii_stub, fext = os.path.splitext(nii_file)
    if fext == '.gz':
        nii_stub, _ = os.path.splitext(nii_stub)

    # Load nifti volume
    print('Opening Nifti-1 volume')
    nii_obj = nib.load(nii_file)

    # Get image header
    hdr = nii_obj.header

    # Image dimensions
    nii_shape = hdr.get_data_shape()
    if len(nii_shape) == 3:
        nx, ny, nz = nii_shape
        nt = 1
    else:
        nx, ny, nz, nt = nii_shape

    # Load image data
    print('Loading voxel data')
    s = nii_obj.get_data()

    print('  Matrix size : (%d, %d, %d, %d)' % (nx, ny, nz, nt))

    # Rescale to 0..255 uint8
    imin, imax = np.min(s), np.max(s)

    print('  Voxel intensity range : [%0.3f, %0.3f]' % (imin, imax))

    s = np.uint8(255.0 * (s - imin) / (imax - imin))

    # Loop over each volume
    for t in range(0, nt):

        # Image stack output directory name
        png_dir = nii_stub + '_%04d' % t

        print('  Volume %d -> %s' % (t,png_dir))

        # Create png_dir if necessary
        if not os.path.exists(png_dir):
            os.makedirs(png_dir)

        # Current volume
        if nt > 1:
            st = s[:,:,:,t]
        else:
            st = s[:,:,:]

        # Loop over each z-slice, outputting as JPEG
        for z in range(0, nz):

            # PNG filename
            png_path = os.path.join(png_dir, png_stub + '_%04d.png' % z)

            # Write single byte image slice to jpg file
            sz_rgb = cv2.cvtColor(st[:, :, z], cv2.COLOR_GRAY2RGB)
            cv2.imwrite(png_path, sz_rgb)

    print('Done')

    # Clean exit
    sys.exit(0)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
