#!/usr/bin/env python3
"""
Output summary statistics for each label in an integer-valued atlas
- volume in voxels and microliters
- center of mass

Author
----
Mike Tyszka, Caltech Brain Imaging Center

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
2019 California Institute of Technology.
"""

__version__ = '0.2.0'

import sys
import argparse
import nibabel as nib
import numpy as np
from scipy.ndimage.measurements import center_of_mass


def main():
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Atlas label statistics')
    parser.add_argument('atlas_file', help="source atlas labels filename")

    args = parser.parse_args()

    atlas_file = args.atlas_file
        
    # Load the source atlas image
    atlas_nii = nib.load(atlas_file)
    atlas_labels = atlas_nii.get_data()
    
    # Atlas voxel volume in mm^3 (microliters)
    atlas_vox_vol_ul = np.array(atlas_nii.header.get_zooms()).prod()
    
    # Create list of unique label values
    labels = np.unique(atlas_labels)

    # Column headers
    print('{:6s} {:10s} {:10s}'.format('Label', 'Voxels', 'ul'))
    
    for label in labels:

        # Skip label 0 (background)
        if label > 0:

            # Extract target label as a boolean mask        
            label_mask = (atlas_labels == label)

            # Center of mass
            com = center_of_mass(label_mask)
            
            # Integrate volume of current label
            label_vol_vox = np.sum(label_mask)
            label_vol_ul = label_vol_vox * atlas_vox_vol_ul
            
            # Only output non-empty labels with index > 0
            if label_vol_ul > 0.0:
                print('%6d %10d %10.1f'.format(label, label_vol_vox, label_vol_ul))


if __name__ == '__main__':
    main()
