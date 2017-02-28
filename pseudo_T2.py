#!/usr/bin/env python3
"""
Generate a placeholder pseduo T2 image from a T1 image
- create signal vs noise mask
- invert signal and mask
- noise conditioning

Authors
----
Mike Tyszka, Caltech Brain Imaging Center

Dates
----
2017-01-12 JMT From scratch

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
2017 California Institute of Technology.
"""

__version__ = '0.1.0'

import sys
import argparse
import numpy as np
import nibabel as nib
from sklearn.cluster import KMeans
from scipy.ndimage.filters import median_filter


def main():

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Simple image segmenter')
    parser.add_argument('-i', '--t1', required=True, help='Input T1 image')

    # Parse command line arguments
    args = parser.parse_args()

    # Get T1w image filename
    t1_fname = args.t1

    # Construct output pseudo T2w filename
    t2_fname = 'pseudo_T2w.nii.gz'

    # Load T1w image
    print('Loading T1w image from %s' % t1_fname)
    t1_nii = nib.load(t1_fname)
    t1_img = t1_nii.get_data()

    # TODO: Write signal masking code
    print('*** Nothing implemented yet ***')

    # Duplicate for now
    print('*** Copying T1w to T2w image ***')
    t2_img = t1_img.copy()

    # Write pseudo T2w image
    print('Saving pseudo T2w image to %s' % t2_fname)
    t2_nii = nib.Nifti1Image(t2_img, t1_nii.get_affine())
    t2_nii.to_filename(t2_fname)

    # Clean exit
    sys.exit(0)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()