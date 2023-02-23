#!/usr/bin/env python3
"""
Merge cortical and subcortical Harvard-Oxford atlas labels into single image
- 1 mm spatial resolution
- Deterministic threshold p > 0.25
- Split cortical labels into LH and RH

AUTHOR
----
Mike Tyszka, Ph.D.

LICENSE
----

MIT License

Copyright (c) 2019 Mike Tyszka

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import os
import sys
import csv
import numpy as np
import nibabel as nib


def main():

    # Flags
    save_intermediates = True

    # FSL Harvard-Oxford source directories and images
    fsl_dir = os.environ['FSL_DIR']
    atlas_dir = os.path.join(fsl_dir, 'data', 'atlases')
    ho_dir = os.path.join(atlas_dir, 'HarvardOxford')
    cort_labels_fname = os.path.join(ho_dir, 'HarvardOxford-cort-maxprob-thr25-1mm.nii.gz')
    subcort_labels_fname = os.path.join(ho_dir, 'HarvardOxford-sub-maxprob-thr25-1mm.nii.gz')

    try:
        cort_nii = nib.load(cort_labels_fname)
        cort_labels = cort_nii.get_data()
    except IOError:
        raise

    try:
        subcort_nii = nib.load(subcort_labels_fname)
        subcort_labels = subcort_nii.get_data()
    except IOError:
        raise

    # Divide cortical labels into left and right hemispheres
    hx = 90
    cort_lh_labels = cort_labels.copy()
    cort_lh_labels[0:hx, :, :] = 0
    cort_rh_labels = cort_labels.copy()
    cort_rh_labels[hx:, :, :] = 0

    # Bilateral cortical labels run from 1 to 48 (0 = background)
    # Add 48 to RH labels
    rh_cort_mask = cort_rh_labels > 0
    cort_rh_labels[rh_cort_mask] += 48

    # Cortical mask for overlap protection
    lh_cort_mask = cort_lh_labels > 0
    cort_mask = lh_cort_mask + rh_cort_mask

    if save_intermediates:

        # Save LH and RH cortical labels
        cort_lh_nii = nib.Nifti1Image(cort_lh_labels, cort_nii.affine)
        fname = os.path.abspath('cort_lh.nii.gz')
        print('Saving LH cortical labels to {}'.format(fname))
        try:
            cort_lh_nii.to_filename(fname)
        except IOError:
            raise

        cort_rh_nii = nib.Nifti1Image(cort_rh_labels, cort_nii.affine)
        fname = os.path.abspath('cort_rh.nii.gz')
        print('Saving RH cortical labels to {}'.format(fname))
        try:
            cort_rh_nii.to_filename(fname)
        except IOError:
            raise

    # Remove GM, WM and CSF labels from sub-cortical data
    mask_lh = (subcort_labels == 1) + (subcort_labels == 2) + (subcort_labels == 3)
    mask_rh = (subcort_labels == 12) + (subcort_labels == 13) + (subcort_labels == 14)
    mask = mask_lh + mask_rh
    subcort_labels[mask] = 0

    # Add 94 to subcortical labels
    subcort_mask = subcort_labels > 0
    subcort_labels[subcort_mask] += 94

    # Remove cortical-subcortical overlap due to 25% deterministic threshold
    subcort_labels[cort_mask] = 0

    if save_intermediates:

        # Save LH and RH cortical labels
        subcort_nii = nib.Nifti1Image(subcort_labels, subcort_nii.affine)
        fname = os.path.abspath('subcort.nii.gz')
        print('Saving sub-cortical labels to {}'.format(fname))
        try:
            subcort_nii.to_filename(fname)
        except IOError:
            raise

    # Combine cortical and subcortical renumbered labels into single integer-valued image
    all_labels = cort_lh_labels + cort_rh_labels + subcort_labels

    # Save combined labels
    all_nii = nib.Nifti1Image(all_labels, cort_nii.affine)
    fname = os.path.abspath('HarvardOxford-all-maxprob-thr25-1mm.nii.gz')
    print('Saving sub-cortical labels to {}'.format(fname))
    try:
        all_nii.to_filename(fname)
    except IOError:
        raise


if '__main__' in __name__:
    main()
