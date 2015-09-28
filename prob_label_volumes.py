#!/opt/local/Library/Frameworks/Python.framework/Versions/2.7/Resources/Python.app/Contents/MacOS/Python
"""
Output volumes of each probabilistic label in an atlas in microliters.
Accepts multiple 4D prob atlas files

Usage
----
prob_label_volumes.py <atlas_file>
prob_label_volumes.py -h

Example
----
>>> prob_label_volumes.py *_probs.nii.gz

Authors
----
Mike Tyszka, Caltech Brain Imaging Center

Dates
----
2015-05-31 JMT Adapt from label_volumes.py

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
import argparse
import nibabel as nib
import numpy as np


def main():
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Probabilistic label volumes in microliters')
    parser.add_argument('prob_files', type=str, nargs='+', help="List of 4D prob label images")

    # Parse command line arguments
    args = parser.parse_args()
    prob_files = args.prob_files
    
    for p_file in prob_files:
        
        # Force absolute path
        p_file = os.path.abspath(p_file)
        
        # Load the source atlas image
        p_nii = nib.load(p_file)
        p = p_nii.get_data()
        nd = p.ndim
        
        # Atlas voxel volume in mm^3 (microliters)
        atlas_vox_vol_ul = np.array(p_nii.header.get_zooms()).prod()
        
        # Treat probabilities as partial volumes and integrate
        if nd == 3:

            V = np.sum(p) * atlas_vox_vol_ul
            print('%0.3f' % V) 

        elif nd == 4:

            nx,ny,nz,nt = p.shape

            for t in range(0,nt):

                V = np.sum(p[:,:,:,t]) * atlas_vox_vol_ul
                print('%0.3f' % V),

            # Final newline
            print          
    
    # Clean exit
    sys.exit(0)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
