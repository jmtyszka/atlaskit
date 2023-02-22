#!/usr/bin/env python3
"""
Separate bilateral prob atlas into hemispheric atlases
- assumes original atlas is centered on the anatomic midline and uses the voxel midplane to split the volumes
- Outputs a new 4D prob atlas with twice as many label volumes

Authors
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
2022 California Institute of Technology.
"""

__version__ = '0.1.0'

import os.path as op
import sys
import argparse
import nibabel as nib
import numpy as np


def main():
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Probabilistic label volumes in microliters')
    parser.add_argument('-i', '--inatlas', required=True, help="Bilateral 4D probabilistic atlas")
    parser.add_argument('-o', '--outatlas', help="Hemispheric 4D probabilistic atlas")

    # Parse command line arguments
    args = parser.parse_args()
    bilat_fname = args.inatlas
    if args.outatlas:
        hemi_fname = op.realpath(args.outatlas)
    else:
        hemi_fname = args.inatlas.replace('.nii.gz', '_LR.nii.gz')

    # Load the bilateral atlas
    print(f'Loading bilateral atlas from {bilat_fname}')
    bilat_nii = nib.load(bilat_fname)
    bilat_img = bilat_nii.get_fdata()

    # Grab 4D image dimensions
    nx, ny, nz, nl = bilat_img.shape

    # Assume data is in RAS order
    # TODO: determine LR anatomic axis from qform (qform_code == 1)
    hx = int(nx/2.0)

    # LH and RH masks
    rh_mask = np.zeros([nx, ny, nz]).astype(int)
    rh_mask[0:hx, ...] = 1
    lh_mask = 1 - rh_mask

    # Make space for hemispheric prob labels
    rh_img = np.zeros_like(bilat_img)
    lh_img = np.zeros_like(bilat_img)

    print(f'Splitting bilateral labels into hemispheric labels')
    for lc in range(nl):
        print(f'  Label {lc}')
        rh_img[..., lc] = bilat_img[..., lc] * rh_mask
        lh_img[..., lc] = bilat_img[..., lc] * lh_mask

    # Concatenate RH and LH prob atlases
    hemi_img = np.concatenate([rh_img, lh_img], axis=3)

    # Create and save new Nifti image
    print(f'Saving hemispheric atlas to {hemi_fname}')
    hemi_nii = nib.Nifti1Image(hemi_img, affine=bilat_nii.affine)
    nib.save(hemi_nii, hemi_fname)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
