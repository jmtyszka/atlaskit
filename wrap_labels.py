# /opt/local/bin/python
"""
Create a skin surrounding multiple manually labeled sections
- Allows convex labels to be defined with fewer sections
- 

Usage
----
skin_labels.py <label image> <index list>
skin_labels.py -h

Example
----
>>> skin_labels.py labels.nii.gz 1 3 4

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
import skimage as sk
import scipy as sp

def main():
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Skin labels')
    parser.add_argument('label_image', help="Labeled volume")
    parser.add_argument('labels', help="Label numbers to skin")

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

            # Create label mask from A and B volumes
            A_mask = (A_labels == label)
            B_mask = (B_labels == label)
    
            # Count voxels in each mask
            nA, nB = np.sum(A_mask), np.sum(B_mask)
            
            # Only calculate stats if labels present in A or B
            if nA > 0 or nB > 0:
                
                # Find intersection and union of A and B masks
                AandB = np.logical_and(A_mask, B_mask)
                AorB = np.logical_or(A_mask, B_mask)
            
                # Count voxels in intersection and union
                nAandB, nAorB = np.sum(AandB), np.sum(AorB)
              
                # Similarity coefficients
                Jaccard = nAandB / float(nAorB)
                Dice = 2.0 * nAandB / float(nA + nB)
                
                # Absolute volumes of label in A and B
                A_vol_ul = np.sum(A_mask) * atlas_vox_vol_ul
                B_vol_ul = np.sum(B_mask) * atlas_vox_vol_ul
    
                print('%8d %8d %8d %10.3f %10.3f %10.3f %10.3f' %
                    (label, nA, nB, A_vol_ul, B_vol_ul, Jaccard, Dice)) 
    
    # Clean exit
    sys.exit(0)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()