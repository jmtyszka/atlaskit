#!/usr/bin/env python3
"""
Combine several labels in a 4D prob atlas by simple addition (probabilistic OR)

Usage
----
prob_or.py -i <4D prob atlas> -o <3D prob atlas> <space separated label list>
prob_or.py -h

Example
----
>>> prob_or.py -i prob_atlas.nii.gz -o merged_labels.nii.gz 1 2 5 10

Authors
----
Mike Tyszka, Caltech Brain Imaging Center

Dates
----
2015-07-29 JMT From scratch
2016-09-27 JMT Clarify argparse help
2023-02-21 JMT Update to python 3.9

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
    parser = argparse.ArgumentParser(description='Probabilistic OR multiple labels from a 4D probabilistic atlas')
    parser.add_argument('-i', '--input', help='Input 4D prob atlas')
    parser.add_argument('-o', '--output', help='Output 3D prob map')
    parser.add_argument('labels', nargs='+', help='Space-separated list of label indices (zero-indexed)')

    # Parse command line arguments
    args = parser.parse_args()

    in_file = args.input if args.input else 'Probabilistic_Atlas.nii.gz'
    out_file = args.output if args.output else 'Probabilistic_OR.nii.gz'

    # Load probabilistic atlas
    print('Loading probabilistic atlas from %s' % in_file)
    in_nii = nib.load(in_file)
    p = in_nii.get_fdata()

    # Probabilistic OR of selected labels
    print('Probabilistic OR of selected labels')
    p_or = np.sum(p[..., args.labels], axis=3)
    
    # Write 4D probabilistic atlas
    print('Saving result to %s' % out_file)
    prob_nii = nib.Nifti1Image(p_or, in_nii.affine)
    prob_nii.to_filename(out_file)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
