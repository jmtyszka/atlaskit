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
import pandas as pd
import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
from skimage.exposure import rescale_intensity
from skimage import color

__version__ = '0.1'


def main():

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Create series of prob label color overlays')
    parser.add_argument('-a', '--atlasdir', required=False, help='Directory containing probabilistic atlas [.]')
    parser.add_argument('-r', '--rois', required=False, help='ROI specification file [<atlasdir>/rois.txt]')
    parser.add_argument('--bilateral', action='store_true', default=False, help='Use bilateral labels')

    # Parse command line arguments
    args = parser.parse_args()

    if args.atlasdir:
        atlas_dir = args.atlasdir
    else:
        atlas_dir = '.'

    if args.rois:
        roi_specfile = args.rois
    else:
        roi_specfile = os.path.join(atlas_dir, 'rois.txt')

    label_keyfile = os.path.join(atlas_dir, 'labels.txt')

    # Use ITK-SNAP label key colors
    atlas_color = True
    
    print('')
    print('----------------------------')
    print('Probabilistic label overlays')
    print('----------------------------')
    print('')

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
    print('Label key         : %s' % label_keyfile)
    print('ROI spec file     : %s' % roi_specfile)
    print('Bilateral         : %s' % ('Yes' if args.bilateral else 'No'))
    print('')

    # Load ROI specs
    if os.path.isfile(roi_specfile):
        rois = load_rois(roi_specfile)
    else:
        print('* ROI specification file %s does not exist' % roi_specfile)
        sys.exit(1)

    # Load label key from atlas directory
    label_key = load_key(label_keyfile)

    # Extract HSV label colors (n_labels x 3 array)
    hsv = label_rgb2hsv(label_key)

    # Load T1w template
    print('< Loading T1w template')
    T1w_fname = os.path.join(atlas_dir, 'T1w_template.nii.gz')
    T1w_nii = nib.load(T1w_fname)
    T1w_img = T1w_nii.get_data()

    # Load T2w template
    print('< Loading T2w template')
    T2w_fname = os.path.join(atlas_dir, 'T2w_template.nii.gz')
    T2w_nii = nib.load(T2w_fname)
    T2w_img = T2w_nii.get_data()

    # Hardwire scaling for CIT168 T1w and T2w templates
    # TODO: generalize or add some command line arguments
    # T1w : 1.0 - 4.0
    # T2w : 0.5 - 4.0
    T1w_img = rescale_intensity(T1w_img, in_range=(1.0, 4.0), out_range=(0.0, 1.0))
    T2w_img = rescale_intensity(T2w_img, in_range=(0.5, 4.0), out_range=(0.0, 1.0))

    # Load the 4D probabilistic atlas
    print('< Loading probabilistic image')
    if args.bilateral:
        prob_fname = os.path.join(atlas_dir, 'prob_atlas_bilateral.nii.gz')
    else:
        prob_fname = os.path.join(atlas_dir, 'prob_atlas.nii.gz')

    p_nii = nib.load(prob_fname)
    p_atlas = p_nii.get_data()

    print('')

    # Count prob labels
    n_labels = p_atlas.shape[3]

    # Loop over all ROIs
    for index, roi in rois.iterrows():

        print('Processing ROI %s' % roi.roiname)

        # Create ROI output directory
        roi_dir = os.path.join(overlay_dir, roi.roiname)
        if not os.path.isdir(roi_dir):
            os.mkdir(roi_dir)

        # Extract ROIs from templates
        T1w_roi = extract_roi(T1w_img, roi)
        T2w_roi = extract_roi(T2w_img, roi)

        # Colorize template ROIs
        T1w_rgb = tint(T1w_roi, hue=0.0, saturation=0.0)
        T2w_rgb = tint(T2w_roi, hue=0.0, saturation=0.0)

        # Initialize the all-label overlay
        all_prob_rgb = np.zeros_like(T1w_rgb)
        all_det_rgb = np.zeros_like(T1w_rgb)

        # Create equivalent montage for all prob labels with varying hues
        for lc in range(0, n_labels):

            # Extract prob label ROI
            prob_roi = extract_roi(p_atlas[:, :, :, lc], roi)

            # Threshold at p > 0.5
            det_roi = (prob_roi > 0.5).astype(int)

            # Hue and saturation for label overlay
            if atlas_color:
                # Pull HSV from ITK-SNAP label key
                hue, sat, val = hsv[lc+1, 0], hsv[lc+1, 1], hsv[lc+1, 2]
            else:
                # Calculate rotating hue
                cyc_freq = 0.99
                hue = float(np.mod(lc * cyc_freq, n_labels)) / n_labels
                sat, val = 1.0, 1.0

            # print('  %d | %0.2f | %0.2f | %0.2f' % (lc, hue, sat, val))

            # Tint the prob label ROI
            prob_rgb = tint(prob_roi, hue=hue, saturation=sat, value=val)
            det_rgb = tint(det_roi, hue=hue, saturation=sat, value=val)

            # Add tinted overlay to running total
            all_prob_rgb += prob_rgb
            all_det_rgb += det_rgb

        # Overlay prob atlas on bg image
        prob_overlay_rgb = composite(all_prob_rgb, T2w_rgb)
        det_overlay_rgb = composite(all_det_rgb, T2w_rgb)

        # Save bg, prob labels and composite images of the ROI
        save_png(roi_dir, 'T1w.png', T1w_rgb)
        save_png(roi_dir, 'T2w.png', T2w_rgb)
        save_png(roi_dir, 'prob.png', all_prob_rgb)
        save_png(roi_dir, 'prob_overlay.png', prob_overlay_rgb)
        save_png(roi_dir, 'det.png', all_det_rgb)
        save_png(roi_dir, 'det_overlay.png', det_overlay_rgb)

    # Clean exit
    sys.exit(0)


def save_png(overlay_dir, png_fname, img_rgb):

    # Create figure and render montage
    fig = plt.figure(figsize=(15, 10), dpi=100)
    plt.imshow(img_rgb, interpolation='none')
    plt.axis('off')
    plt.legend()

    # Save figure to PNG
    png_path = os.path.join(overlay_dir, png_fname)
    plt.savefig(png_path, bbox_inches='tight')

    # Close figure
    plt.close()


def load_rois(roi_specfile):
    """
    Load ROI specs from specfile

    Header row: roiname x y z dx dy dz
    Row format: %s %d %d %d %d %d %d

    x,y,z is the ROI corner closest to the origin
    dx,dy,dz are the ROI dimensions

    Each ROI should have exactly one of dx,dy or dz = 1 (slice)

    Parameters
    ----------
    roi_specfile

    Returns
    -------

    """

    with open(roi_specfile) as fd:
        rois = pd.read_table(roi_specfile, sep=" ")

    return rois


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


def extract_roi(img, roi, flip_x=False, flip_y=True, flip_z=True):
    """
    Create a montage of all coronal (XZ) slices from a 3D image

    Parameters
    ----------
    img: numpy array
        3D image to montage
    roi: pandas table row
    flip_x: bool
    flip_y: bool
    flip_z: bool

    Returns
    -------

    img_roi: numpy array
    """

    # Source image dimensions
    nx, ny, nz = img.shape

    # Create slices for each axis
    xx = slice(roi.x, roi.x + roi.dx)
    yy = slice(roi.y, roi.y + roi.dy)
    zz = slice(roi.z, roi.z + roi.dz)

    # Extract ROI from volume
    img_roi = img[xx,yy,zz]

    # Squeeze out singlet dimensions
    img_roi = img_roi.squeeze()

    # Construct montage of coronal sections
    return img_roi


def tint(image, hue=0.0, saturation=1.0, value=1.0):
    """
    Add color of the given hue to an RGB image

    Parameters
    ----------
    image: numpy array
    hue: float
    saturation: float
    value: float

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

    # Clamp values [0..1)
    composite_rgb = composite_rgb.clip(0.0, 1.0)

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
