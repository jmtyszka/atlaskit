# /opt/local/bin/python
"""
Interpolate label between sparse sections
- Speeds up manual labeling for larger label volumes

Usage
----
interp_labels.py <label image> <index list>
interp_labels.py -h

Example
----
>>> interp_labels.py labels.nii.gz 1 3 4

Authors
----
Mike Tyszka, Caltech Brain Imaging Center

Dates
----
2015-09-28 JMT From scratch

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
import nibabel as nib
import numpy as np
from scipy.interpolate import Rbf

def main():
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Interpolate labels')
    parser.add_argument('label_image', help="Labeled volume")
    parser.add_argument('labels', help="Label numbers to interpolate")

    # Parse command line arguments
    args = parser.parse_args()
    label_fname = args.label_image

    # Load labeled volume
    label_nii = nib.load(label_fname)
    labels = label_nii.get_data()
    
    if args.labels:
        label_nos = args.labels
    else:
        # Construct list of unique label values in image
        label_nos = np.unique(labels)

    # loop over each unique label value
    for label in label_nos:
        
        if label > 0:

            # Extract current label
            L = int(labels == label)
            
            # Find coordinates of all True voxels
            x, y, z = np.where(L)
            
            # Count number of voxels in label
            nvox = np.sum(L)
            
            # All node values are one
            d = np.ones([nvox])
            
            # RBF interpolate
            rbf = Rbf(x,y,z,d)
    
            # Interpolated label
            Li = rbf(xi,yi,zi)
    
    # Clean exit
    sys.exit(0)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()