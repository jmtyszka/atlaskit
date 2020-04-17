#!/usr/bin/env python3
"""
Calculate an optionally padded ROI for a single probabilistic label

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
2020 California Institute of Technology.
"""

__version__ = '0.1.0'

import os
import sys
import argparse
import nibabel as nib
import numpy as np


def main():
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Padded ROI for a single probabilistic label')
    parser.add_argument('-p', '--problabel',
                        required=True,
                        help='Probabilistic label file')
    parser.add_argument('-i', '--infiles',
                        required=False,
                        default=[],
                        action='append',
                        help='Image to crop to ROI')
    parser.add_argument('-v', '--voxelpad',
                        required=False,
                        default='8',
                        help='Voxel padding [8]')
    parser.add_argument('-t', '--threshold',
                        required=False,
                        default='0.1',
                        help='Probability threshold [0.1]')

    # Parse command line arguments
    args = parser.parse_args()

    # Collect arguments
    plabel_fname = os.path.abspath(args.problabel)
    vox_pad = np.int(args.voxelpad)
    prob_thresh = np.float(args.threshold)

    # Load the prob label image
    p_nii = nib.load(plabel_fname)
    p = p_nii.get_data()

    # Original affine tx and voxel dimensions
    p_affine = p_nii.affine
    p_hdr = p_nii.header
    p_vx, p_vy, p_vz = p_hdr.get_zooms()

    # Threshold at p > threshold
    mask = p > prob_thresh

    # Find minimum bounding box at this threshold
    xmin, xmax, ymin, ymax, zmin, zmax = bounding_box(mask)

    # Apply padding
    xmin, xmax = xmin - vox_pad, xmax + vox_pad
    ymin, ymax = ymin - vox_pad, ymax + vox_pad
    zmin, zmax = zmin - vox_pad, zmax + vox_pad

    # Clamp range to image volume
    ny, nx, nz = mask.shape
    xmin = np.clip(xmin, 0, nx)
    xmax = np.clip(xmax, 0, nx)
    ymin = np.clip(ymin, 0, ny)
    ymax = np.clip(ymax, 0, ny)
    zmin = np.clip(zmin, 0, nz)
    zmax = np.clip(zmax, 0, nz)

    # Bounding box dimensions
    dx = xmax - xmin + 1
    dy = ymax - ymin + 1
    dz = zmax - zmin + 1

    # Update transform matrix to account for crop
    crop_affine = p_affine.copy()
    crop_affine[0:3, 3] += np.array([xmin * p_vx, ymin * p_vy, zmin * p_vz])

    # Crop any images provided to the bounding box
    for img_fname in args.infiles:

        img_nii = nib.load(img_fname)
        img = img_nii.get_data()
        ndims = np.ndim(img)

        if ndims == 3:
            img_mask = img[xmin:xmax, ymin:ymax, zmin:zmax]
        elif ndims == 4:
            img_mask = img[xmin:xmax, ymin:ymax, zmin:zmax, :]
        else:
            print('* Unsupported image dimensionality {}'.format(ndims))
            sys.exit(1)

        img_mask_nii = nib.Nifti1Image(img_mask, affine=crop_affine)

        img_mask_outfile = img_fname.replace('.nii.gz', '_crop.nii.gz')
        img_mask_nii.to_filename(img_mask_outfile)

    # Output the bounding box limits in fslmaths/fslroi format
    # xmin dx ymin dy zmin dz tmin dt
    print(xmin, dx, ymin, dy, zmin, dz, 0, 1)


def bounding_box(img):

    # Collapse in z and y - sufficient for bounding box
    xyproj = np.sum(img, axis=2)
    xzproj = np.sum(img, axis=1)

    # Collapse 2D projections in x, y and z
    xproj = np.sum(xyproj, axis=1)
    yproj = np.sum(xyproj, axis=0)
    zproj = np.sum(xzproj, axis=0)

    xmin, xmax = np.where(xproj)[0][[0, -1]]
    ymin, ymax = np.where(yproj)[0][[0, -1]]
    zmin, zmax = np.where(zproj)[0][[0, -1]]

    return xmin, xmax, ymin, ymax, zmin, zmax


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
