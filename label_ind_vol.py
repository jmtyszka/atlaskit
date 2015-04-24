#!/opt/local/Library/Frameworks/Python.framework/Versions/2.7/Resources/Python.app/Contents/MacOS/Python
"""
Estimate label volumes for an individual using the atlas label volume and
local jacobian determinant of the template to individual mapping.

Usage
----
label_ind_vol.py <atlas> <individual jacobian map>
label_ind_vol.py -h

Example
----
>>> label_ind_vol.py atlas.nii.gz 105014_temp2ind_detjac.nii.gz

Authors
----
Mike Tyszka, Caltech Brain Imaging Center

Dates
----
2015-04-13 JMT From scratch

License
----
This file is part of mrgaze.

    atlaskit is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    atlaskit is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with mrgaze.  If not, see <http://www.gnu.org/licenses/>.

Copyright
----
2015 California Institute of Technology.
"""

__version__ = '0.1.0'

import sys
import argparse
import nibabel as nib
import numpy as np


def main():
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Individual atlas label volumes')
    parser.add_argument('atlas_file', help="source atlas labels filename")
    parser.add_argument('jacobian_file', help="local jacobian determinant map")
    parser.add_argument('csv_file', help="individual label volumes CSV file")

    args = parser.parse_args()

    atlas_file = args.atlas_file
    jac_file = args.jacobian_file
    csv_file = args.csv_file
        
    # Load the source atlas image
    print('Loading atlas labels from %s' % atlas_file)
    atlas_nii = nib.load(atlas_file)
    atlas_labels = atlas_nii.get_data()
    
    # Atlas voxel volum in mm^3 (microliters)
    atlas_vox_vol_ul = np.array(atlas_nii.header.get_zooms()).prod()
    
    print('Atlas voxel volume : %0.5f ul' % atlas_vox_vol_ul)
    
    # Create list of label values (not necessarily contiguous)
    max_l = atlas_labels.max()
    labels = range(0, max_l+1)        
    
    # Load jacobian image
    print('Loading Jacobian image from %s' % jac_file)
    jac_nii = nib.load(jac_file)
    jac_img = jac_nii.get_data()
    
    # Make space for label volumes
    label_vol_ul = np.zeros([max_l+1], dtype=int)

    for label in labels:

        # Extract target label as a boolean mask        
        label_mask = (atlas_labels == label)

        # Integrate local Jacobian determinant over mask
        total_jac = jac_img[label_mask].sum()
        label_vol_ul[label] = total_jac * atlas_vox_vol_ul
        
        print('  Label %03d  : %0.1f ul' % (label, label_vol_ul[label])) 
    
    # Save label volumes to CSV text file
    print('Saving label volumes to %s' % csv_file)
    np.savetxt(csv_file, label_vol_ul, delimiter=',', fmt='%0.1f')
    
    print('Done')
    
    # Clean exit
    sys.exit(0)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()