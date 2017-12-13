#!/usr/bin/env python
"""
Create a series of color overlays of prob labels on cor, sag or ax ROIs
- ROIs are specified in a text file

Usage
----
atlas_report.py -a <atlas directory created by atlas.py>
atlas_report.py -h

Authors
----
Mike Tyszka, Caltech Brain Imaging Center

Dates
----
2017-02-21 JMT Split from atlas.py

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
import jinja2
import colorsys
import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
from datetime import datetime
from skimage.util.montage import montage2d
from skimage import color
from atlas import get_label_name

__version__ = '0.1'


def main():

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Create series of prob label color overlays')
    parser.add_argument('-a', '--atlasdir', required=True, help='Directory containing probabilistic atlas')
    parser.add_argument('-r', '--rois', required=True, help='ROI specification file')

    # Parse command line arguments
    args = parser.parse_args()
    atlas_dir = args.atlasdir
    roi_specfile = args.rois

    # Use ITK-SNAP label key colors
    atlas_color = False
    
    print('')
    print('----------------------------')
    print('Probabilistic label overlays')
    print('----------------------------')

    # Check for atlas directory existence
    if not os.path.isdir(atlas_dir):
        print('Atlas directory does not exist (%s) - exiting' % atlas_dir)
        sys.exit(1)

    # Create overlay directory within atlas directory
    overlay_dir = os.path.join(atlas_dir, 'overlays')
    if not os.path.isdir(overlay_dir):
        os.mkdir(overlay_dir)

    print('Atlas directory   : %s' % atlas_dir)
    print('Overlay directory : %s' % overlay_dir)
    print('')

    # Load ROI specs
    if os.path.isfile(roi_specfile):
        rois = load_rois(roi_specfile)
    else:
        print('* ROI specification file %s does not exist' % roi_specfile)
        sys.exit(1)

    # CIT atlas directory from shell environment
    cit_dir = os.environ['CIT168_DIR']

    # Load label key from atlas directory
    label_key = load_key(os.path.join(atlas_dir, 'labels.txt'))

    # Extract HSV label colors (n_labels x 3 array)
    hsv = label_rgb2hsv(label_key)

    # Load background image
    print('  Loading background image')
    bg_fname = os.path.join(cit_dir, 'CIT168_700um', 'CIT168_T1w_700um.nii.gz')
    bg_nii = nib.load(bg_fname)
    bg_img = bg_nii.get_data()

    # Normalize background intensity range to [0,1]
    bg_img = bg_img / np.max(bg_img)

    # Load the 4D probabilistic atlas
    print('  Loading probabilistic image')
    p_nii = nib.load(os.path.join(atlas_dir, overlay_fname))
    p_atlas = p_nii.get_data()

    # Count prob labels
    n_labels = p_atlas.shape[3]

    # Loop over all rois
    for roi in rois:

        # Create montage of coronal sections through cropped bg image
        bg_mont = coronal_montage(bg_crop, n_rows, n_cols)

        bg_mont_rgb = tint(bg_mont, hue=0.0, saturation=0.0)

        # Initialize the all-label overlay
        overlay_mont_rgb = np.zeros_like(bg_mont_rgb)

        # Create equivalent montage for all prob labels with varying hues
        for lc in range(0, n_labels):

            # Construct prob label montage
            p_mont = coronal_montage(p_crop[:, :, :, lc], n_rows, n_cols)

            # Hue and saturation for label overlay
            if atlas_color:
                # Pull HSV from ITK-SNAP label key
                hue, sat, val = hsv[lc, 0], hsv[lc, 1], hsv[lc, 2]
            else:
                # Calculate rotating hue
                hue = float(np.mod(lc * 3, n_labels)) / n_labels
                sat, val = 1.0, 1.0

            # Tint the montage
            p_mont_rgb = tint(p_mont, hue=hue, saturation=sat, value=val)

            # Add tinted overlay to running total
            overlay_mont_rgb += p_mont_rgb

        # Composite prob atlas overlay on bg image
        mont_rgb = composite(overlay_mont_rgb, bg_mont_rgb)

        # Create figure and render montage
        fig = plt.figure(figsize=(15, 10), dpi=100)
        plt.imshow(mont_rgb, interpolation='none')
        plt.axis('off')
        plt.legend()

        # Save figure to PNG
        montage_fname = overlay_fname.replace('.nii.gz', '_montage.png')
        print('  Saving image to %s' % montage_fname)
        plt.savefig(os.path.join(report_dir, montage_fname), bbox_inches='tight')

    # Clean exit
    sys.exit(0)


def label_rgb2hsv(label_key):
    """
    Extract label RGB colors and convert to HSV

    Parameters
    ----------
    label_key: data frame

    Returns
    -------
    hsv: numpy array

    """

    rgb = np.array(label_key[['R','G','B']]) / 255.0
    rgb = rgb.reshape([rgb.shape[0], 1, 3])
    hsv = color.rgb2hsv(rgb)
    hsv = hsv.reshape([-1,3])

    return hsv


def coronal_montage(img, n_rows=4, n_cols=4, flip_x=False, flip_y=True, flip_z=True):
    """
    Create a montage of all coronal (XZ) slices from a 3D image

    Parameters
    ----------
    img: 3D image to montage
    n_rows: number of montage rows
    n_cols: number of montage columns
    rot: CCW 90deg rotations to apply to each section

    Returns
    -------

    cor_mont: coronal slice montage of img
    """

    # Total number of sections to extract
    n = n_rows * n_cols

    # Source image dimensions
    nx, ny, nz = img.shape

    # Coronal (XZ) sections
    yy = np.linspace(0, ny-1, n).astype(int)
    cors = img[:,yy,:]

    if flip_x:
        cors = np.flip(cors, axis=0)
    if flip_y:
        cors = np.flip(cors, axis=1)
    if flip_z:
        cors = np.flip(cors, axis=2)

    # Permute image axes for montage2d: original y becomes new x
    img = np.transpose(cors, (1,2,0))

    # Construct montage of coronal sections
    cor_mont = montage2d(img, fill=0, grid_shape=(n_rows, n_cols))

    return cor_mont


def tint(image, hue=0.0, saturation=1.0, value=1.0):
    """
    Add color of the given hue to an RGB image

    Parameters
    ----------
    image: numpy array
    hue: float
    saturation: float

    Returns
    -------
    """

    hsv = np.zeros([image.shape[0], image.shape[1], 3])
    hsv[:, :, 0] = hue
    hsv[:, :, 1] = saturation
    hsv[:, :, 2] = image * value

    return color.hsv2rgb(hsv)


def composite(overlay_rgb, background_rgb):
    """
    Alpha composite RGB overlay on RGB background
    - derive alpha from HSV value of overlay

    Parameters
    ----------
    overlay_rgb:
    background_rgb:

    Returns
    -------

    """

    overlay_hsv = color.rgb2hsv(overlay_rgb)
    value = overlay_hsv[:,:,2]
    alpha_rgb = np.dstack((value, value, value))

    composite_rgb = overlay_rgb * alpha_rgb + background_rgb * (1.0 - alpha_rgb)

    return composite_rgb


def load_key(key_fname):
    """
    Parse an ITK-SNAP label key file

    Parameters
    ----------
    key_fname: ITK-SNAP label key filename

    Returns
    -------
    key: Data table containing ITK-SNAP style label key
    """

    import pandas as pd

    # Import key as a data table
    # Note the partially undocumented delim_whitespace flag
    key = pd.read_table(key_fname,
                        comment='#',
                        header=None,
                        names=['Index', 'R', 'G', 'B', 'A', 'Vis', 'Mesh', 'Name'],
                        delim_whitespace=True)

    return key


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    try:
        cit_dir = os.environ['CIT168_DIR']
    except KeyError:
        print('* Environmental variable CIT168_DIR not set - exiting')
        sys.exit(1)
        
    main()
