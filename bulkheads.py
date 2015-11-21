#!/opt/local/bin/python
'''
Find bulkhead edge voxels in sparsely labeled volume

Authors
----
Mike Tyszka, Caltech Brain Imaging Center

Dates
----
2015-11-10 JMT From scratch

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
'''

__version__ = '0.1.0'

import os
import sys
import argparse
import nibabel as nib
import numpy as np
from scipy.ndimage.filters import convolve


def main():
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Interpolate labels')
    parser.add_argument('-i','--input', required=True, help="Labeled volume")

    # Parse command line arguments
    args = parser.parse_args()

    # Get mandatory filename argument
    label_fname = args.input

    # Construct output filename    
    out_stub, out_ext = os.path.splitext(label_fname)
    if out_ext == '.gz':
        out_stub, _ = os.path.splitext(out_stub)
    out_fname = out_stub + '_bulkheads.nii.gz'

    # Load labeled volume
    label_nii = nib.load(label_fname)
    labels = label_nii.get_data()
    
    # Size of image space
    nx, ny, nz = labels.shape
    
    # Destination label volume
    bulkhead = np.zeros_like(labels, dtype='float')
    
    # Construct list of unique label values in image
    label_nos = np.unique(labels)
    
    # Convolution filter weights
    w = np.ones([3,3,3]) * -1.0
    w[1,1,1] = 1.0
    
    # loop over each unique label value
    for label in label_nos:
        
        if label > 0:
            
            print('Find bulkheads in label %d' % label)

            # Extract current label
            L = (labels == label).astype(float)
            
            # 3D Sobel filter label image
            bulkhead = convolve(L, w)
            
    
    # Save interpolated label volume
    print('Saving bulkhead edge image to %s' % out_fname)
    out_nii = nib.Nifti1Image(bulkhead, label_nii.get_affine())
    out_nii.to_filename(out_fname)
        
    
    # Clean exit
    sys.exit(0)

    

    

# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
