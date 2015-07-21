#!/opt/local/Library/Frameworks/Python.framework/Versions/2.7/Resources/Python.app/Contents/MacOS/Python
"""
Calculate Dice, Jaquard and related stats for inter and intra-observer labeling comparisons

Usage
----
dice.py <labels_A> <labels_B>
dice.py -h

Example
----
>>> dice.py labels_A.nii.gz labels_B.nii.gz

Authors
----
Mike Tyszka, Caltech Brain Imaging Center

Dates
----
2015-07-21 JMT From scratch

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


def main():
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Similarity metrics for labeled volumes')
    parser.add_argument('labels_A', help="Labeled volume A")
    parser.add_argument('labels_B', help="Labeled volume B")

    # Parse command line arguments
    args = parser.parse_args()
    labels_A = args.labels_A
    labels_B = args.labels_B
        
    # Load labeled volumes
    A_nii, B_nii = nib.load(labels_A), nib.load(labels_B)
    A_labels, B_labels = A_nii.get_data(), B_nii.get_data()
    
    # Voxel volume in mm^3 (microliters) (from volume A)
    atlas_vox_vol_ul = np.array(A_nii.header.get_zooms()).prod()
    
    # Create list of possible label values from volume A
    max_l = np.int(A_labels.max())
    labels = range(0, max_l+1)

    # Colume headers
    print('%8s %8s %8s %10s %10s %10s %10s' %
        ('Label', 'nA', 'nB', 'vA_ul', 'vB_ul', 'Jaccard', 'Dice')) 
    
    for label in labels:

        # Skip label 0 (background)
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
