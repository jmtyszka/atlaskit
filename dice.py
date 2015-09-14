#!/opt/local/bin/python
"""
Calculate Dice, Jaquard and related stats for inter and intra-observer labeling comparisons

Usage
----
dice.py <labelsA> <labelsB>
dice.py -h

Example
----
>>> dice.py labelsA.nii.gz labelsB.nii.gz

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

__version__ = '0.2.0'

import sys
import argparse
import nibabel as nib
import numpy as np


def main():
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Similarity metrics for labeled volumes')
    parser.add_argument('labelsA', help="Labeled volume A")
    parser.add_argument('labelsB', help="Labeled volume B")


    # Parse command line arguments
    args = parser.parse_args()
    labelsA, labelsB = args.labelsA, args.labelsB
        
    # Load labeled volumes
    A_nii, B_nii = nib.load(labelsA), nib.load(labelsB)
    A_labels, B_labels = A_nii.get_data(), B_nii.get_data()
    
    # Voxel volume in mm^3 (microliters) (from volume A)
    atlas_vox_vol_ul = np.array(A_nii.header.get_zooms()).prod()

    # Colume headers
    print('%8s,%8s,%8s,%10s,%10s,%10s,%10s,' %
        ('Label', 'nA', 'nB', 'vA_ul', 'vB_ul', 'Jaccard', 'Dice')) 
    
    # Construct list of unique label values in image
    unique_labels = np.unique(A_labels)

    # loop over each unique label value
    for label in unique_labels:
        
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
    
                print('%8d,%8d,%8d,%10.3f,%10.3f,%10.3f,%10.3f,' %
                    (label, nA, nB, A_vol_ul, B_vol_ul, Jaccard, Dice)) 
    
    # Clean exit
    sys.exit(0)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
