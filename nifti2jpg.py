#!/usr/bin/env python3
"""
Convert single Nifti-1 3D or 4D volume to set of JPEG stacks.

Usage
----
nifti2jpg.py -i <Nifti filename> -o <JPEG image stub>
nifti2jpg.py -h

Example
----
>>> nifti2jpg.py -i mynifti.nii.gz -o myjpgs 

Authors
----
Mike Tyszka, Caltech Brain Imaging Center

Dates
----
2015-09-03 JMT From scratch

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

import os
import sys
import cv2
import argparse
import nibabel as nib
import numpy as np

def main():
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Convert 3D Nifti volume to JPEG stack')
    parser.add_argument('-i','--nii_file', help='3D or 4D Nifti image')
    parser.add_argument('-o','--jpg_stub', help='JPEG stack stub')
    
    args = parser.parse_args()
    
    nii_file = args.nii_file
    jpg_stub = args.jpg_stub
    
    # Strip extension from Nifti filename for use as a stub
    nii_stub, fext = os.path.splitext(nii_file)
    if fext == '.gz':
        nii_stub, _ = os.path.splitext(nii_stub)
    
    # Load nifti volume
    print('Opening Nifti-1 volume')
    nii_obj = nib.load(nii_file)
            
    # Image dimensions
    nx,ny,nz,nt = nii_obj.header.get_data_shape()

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
        jpg_dir = nii_stub + '_%04d' % t

        print('  Volume %d -> %s' % (t, jpg_dir))
        
        # Create jpg_dir if necessary
        if not os.path.exists(jpg_dir):
            os.makedirs(jpg_dir)
            
        # Current volume
        st = s[:,:,:,t]
        
        # Loop over each z-slice, outputting as JPEG
        for z in range(0, nz):
    
            # JPEG filename
            jpg_path = os.path.join(jpg_dir, jpg_stub + '_%04d.jpg' % z)
            
            # Write single byte image slice to jpg file
            sz_rgb = cv2.cvtColor(st[:,:,z], cv2.COLOR_GRAY2RGB)
            cv2.imwrite(jpg_path, sz_rgb)

    print('Done')
    
    # Clean exit
    sys.exit(0)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
