#!/usr/bin/env python3
"""
Analyse absolute and relative volumes of probabilistic atlas labels occupied by one or
more lesion labels

Usage
----
lesion_atlas_analysis.py
    -l <3D lesion label image>
    -a <4D prob atlas image>
    [-lk <lesion label key>]
    [-ak <atlas label key>]

Authors
----
Mike Tyszka, Caltech Brain Imaging Center

Dates
----
2016-10-26 JMT From scratch

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
2016 California Institute of Technology.
"""

__version__ = '0.1.0'

import os
import sys
import argparse
import nibabel as nib
import numpy as np
from nipype.interfaces import ants


def main():

    # Construct a command line argument parser
    parser = argparse.ArgumentParser(description='Atlas-based lesion volumetrics')
    parser.add_argument('-l', '--lesion', required=True, help='3D multivalued lesion label image')
    parser.add_argument('-a', '--atlas', required=True, help='4D probabilistic atlas label image')
    parser.add_argument('-lk', '--lesionkey', required=False, help='Lesion label key (ITKSNAP format)')
    parser.add_argument('-ak', '--atlaskey', required=False, help='Atlas label key (ITKSNAP format)')

    # Parse command line arguments
    args = parser.parse_args()
    lesion_fname = args.lesion
    atlas_fname = args.atlas

    if args.lesionkey:
        lesion_key = args.lesionkey
    else:
        lesion_key = []

    if args.atlaskey:
        atlas_key = args.lesionkey
    else:
        atlas_key = []

    # Load lesion label image
    try:
        lesion_nii = nib.load(lesion_fname)
        lesion = lesion_nii.
    except:
        print('* Problem loading lesion image')

    # Load probabilistic atlas
    try:
        atlas_nii = nib.load(atlas_fname)
    except:
        print('* Problem loading atlas image')




    # Clean exit
    sys.exit(0)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
