#!/usr/bin/env python3
"""
Create GM, WM and CSF tissue masks from Freesurfer output
 - requires Freesurfer reconstructed subject directory
 - uses mri/ribbon.mgz, wmparc.mgz and T1.mgz

Authors
----
Mike Tyszka, Caltech Brain Imaging Center

Dates
----
2017-06-21 JMT Expand from Julien Dubois code fragment

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


import os
import sys
import argparse
import nibabel as nib
import numpy as np

def main():

    # Cerebellum flag
    include_cerebellum = False

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Create GM, WM and CSF masks from Freesurfer parcellation')
    parser.add_argument('-s', '--subjid', required=True, help='Freesurfer subject ID')
    parser.add_argument('-o', '--outdir', default='.', help='Output directory for tissue masks')

    args = parser.parse_args()
    subjid = args.subjid
    out_dir = args.outdir

    # Freesurfer subject recon directory
    # Get Freesurfer subjects directory from environment ($SUBJECTS_DIR)
    subj_dir = os.path.join(os.environ['SUBJECTS_DIR'], subjid)

    ribbon_fname = os.path.join(subj_dir, 'mri', 'ribbon.mgz')
    wmparc_fname = os.path.join(subj_dir, 'mri', 'wmparc.mgz')
    t1_fname = os.path.join(subj_dir, 'mri', 'T1.mgz')

    # Load FS parcellations
    try:
        print('+ Loading GM ribbon')
        ribbon_mgz = nib.load(ribbon_fname)
        ribbon_img = ribbon_mgz.get_data()
    except:
        print('* Problem opening %s' % ribbon_fname)
        sys.exit(1)

    try:
        print('+ Loading WM parcellation')
        wmparc_mgz = nib.load(wmparc_fname)
        wmparc_img = wmparc_mgz.get_data()
    except:
        print('* Problem opening %s' % wmparc_fname)
        sys.exit(1)

    # Load FS T1 reference
    try:
        print('+ Loading reference T1')
        t1_mgz = nib.load(t1_fname)
        t1_img = t1_mgz.get_data()
    except:
        print('* Problem opening %s' % t1_fname)
        sys.exit(1)

    # Parcel lists for masks

    # Left-Cerebral-White-Matter, Right-Cerebral-White-Matter
    ribbonWMstructures = [2, 41]

    # Left-Cerebral-Cortex, Right-Cerebral-Cortex
    ribbonGMstructures = [3, 42]

    # Fornix, CC-Posterior, CC-Mid-Posterior, CC-Central, CC-Mid-Anterior, CC-Anterior
    wmparcCCstructures = [250, 251, 252, 253, 254, 255]

    # Left-Lateral-Ventricle, Left-Inf-Lat-Vent, 3rd-Ventricle, 4th-Ventricle, CSF
    # Left-Choroid-Plexus, Right-Lateral-Ventricle, Right-Inf-Lat-Vent, Right-Choroid-Plexus
    wmparcCSFstructures = [4, 5, 14, 15, 24, 31, 43, 44, 63]

    if include_cerebellum:

        # Cerebellar-White-Matter-Left, Brain-Stem, Cerebellar-White-Matter-Right
        wmparcWMstructures = [7, 16, 46]

        # Left-Cerebellar-Cortex, Right-Cerebellar-Cortex, Thalamus-Left, Caudate-Left
        # Putamen-Left, Pallidum-Left, Hippocampus-Left, Amygdala-Left, Accumbens-Left
        # Diencephalon-Ventral-Left, Thalamus-Right, Caudate-Right, Putamen-Right
        # Pallidum-Right, Hippocampus-Right, Amygdala-Right, Accumbens-Right
        # Diencephalon-Ventral-Right
        wmparcGMstructures = [8, 47, 10, 11, 12, 13, 17, 18, 26, 28, 49, 50, 51, 52, 53, 54, 58, 60]

    else:

        # Omit cerebellum and brain stem
        wmparcWMstructures = []
        wmparcGMstructures = [10, 11, 12, 13, 17, 18, 26, 28, 49, 50, 51, 52, 53, 54, 58, 60]

    # Construct masks
    wm_mask = np.double(
        np.logical_and(
            np.logical_and(
                np.logical_or(
                    np.logical_or(
                        np.in1d(ribbon_img, ribbonWMstructures), np.in1d(wmparc_img, wmparcWMstructures)),
                    np.in1d(wmparc_img, wmparcCCstructures)),
                np.logical_not(np.in1d(wmparc_img, wmparcCSFstructures))),
            np.logical_not(np.in1d(wmparc_img, wmparcGMstructures))))

    csf_mask = np.double(np.in1d(wmparc_img, wmparcCSFstructures))
    gm_mask = np.double(np.logical_or(np.in1d(ribbon_img,ribbonGMstructures), np.in1d(wmparc_img,wmparcGMstructures)))

    # Reshape mask to 3D
    wm_mask = np.reshape(wm_mask, ribbon_img.shape)
    csf_mask = np.reshape(csf_mask, ribbon_img.shape)
    gm_mask = np.reshape(gm_mask, ribbon_img.shape)

    # Save masks
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    wm_fname = os.path.join(out_dir, 'fs_wm.nii.gz')
    print('+ Saving WM mask to %s' % wm_fname)
    prob_nii = nib.Nifti1Image(wm_mask, ribbon_mgz.get_affine())
    prob_nii.to_filename(wm_fname)

    csf_fname = os.path.join(out_dir, 'fs_csf.nii.gz')
    print('+ Saving CSF mask to %s' % csf_fname)
    prob_nii = nib.Nifti1Image(csf_mask, ribbon_mgz.get_affine())
    prob_nii.to_filename(csf_fname)

    gm_fname = os.path.join(out_dir, 'fs_gm.nii.gz')
    print('+ Saving GM mask to %s' % gm_fname)
    prob_nii = nib.Nifti1Image(gm_mask, ribbon_mgz.get_affine())
    prob_nii.to_filename(gm_fname)

    t1_fname = os.path.join(out_dir, 'fs_t1.nii.gz')
    print('Saving reference T1 to %s' % t1_fname)
    prob_nii = nib.Nifti1Image(t1_img, ribbon_mgz.get_affine())
    prob_nii.to_filename(t1_fname)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()