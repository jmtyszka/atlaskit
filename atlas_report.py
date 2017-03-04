#!/usr/bin/env python3
"""
Create a report of intra and inter-observer atlas label statistics
- requires that atlas.py has been run previously on the labels directory

Usage
----
atlas_report.py -d <labels directory> [-k <ITK-SNAP label key>]
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

__version__ = '0.1.0'

import os
import sys
import csv
import argparse
import jinja2
import numpy as np
import pandas as pd
import nibabel as nib
import matplotlib.pyplot as plt
from glob import glob
from scipy.ndimage import measurements, find_objects
from skimage.io import imsave

def main():

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Create labeling report for a probabilistic atlas')
    parser.add_argument('-a','--atlasdir', required=True, help='Directory containing probabilistic atlas')

    # Parse command line arguments
    args = parser.parse_args()
    atlas_dir = args.atlasdir
    print('Atlas directory : %s' % atlas_dir)

    # Check for atlas directory existence
    if not os.path.isdir(atlas_dir):
        print('Atlas directory does not exist (%s) - exiting' % atlas_dir)
        sys.exit(1)

    # Create report directory within atlas directory
    report_dir = os.path.join(atlas_dir, 'report')
    if not os.path.isdir(report_dir):
        os.mkdir(report_dir)
    print('Report directory : %s' % report_dir)

    print('Loading similarity metrics')
    intra_stats, inter_stats = load_metrics(atlas_dir)

    # Intra-observer reports (one per observer)
    print('Generating intra-observer report')
    intra_observer_reports(report_dir, intra_stats)

    # Inter-observer report
    print('Generating inter-observer report')
    inter_observer_report(report_dir, inter_stats)

    # Summary report page
    print('Writing report summary page')
    summary_report(atlas_dir, intra_stats)

    # Clean exit
    sys.exit(0)


def summary_report(atlas_dir, intra_stats):
    """
    Summary report for the entire atlas
    - projections with overlays

    Parameters
    ----------
    atlas_dir: report directory path
    intra_stats: numpy array of intra-observer stats

    Returns
    -------

    """

    report_dir = os.path.join(atlas_dir, 'report')

    # Intra stats columns: idx, trueidx, obs, tmpA, tmpB, dice, hausdorff
    n_obs = len(np.unique(intra_stats[:,2]))
    n_tmp = len(np.unique(intra_stats[:,3]))

    # Create tryptic overlays through CoM of label
    tryptic_overlays(atlas_dir)

    # Setup template
    template_loader = jinja2.FileSystemLoader(searchpath="/Users/jmt/GitHub/atlaskit")
    template_env = jinja2.Environment(loader=template_loader)
    template_fname = "atlas_summary.jinja"
    template = template_env.get_template(template_fname)

    # Template variables
    template_vars = {"n_obs": n_obs,
                     "n_tmp": n_tmp}

    # Finally, process the template to produce our final text.
    output_text = template.render(template_vars)

    # Write page to report directory
    with open(os.path.join(os.path.join(atlas_dir), 'report', 'index.html'), "w") as f:
        f.write(output_text)


def tryptic_overlays(atlas_dir):
    """

    Parameters
    ----------
    atlas_dir

    Returns
    -------

    """

    report_dir = os.path.join(atlas_dir, 'report')

    # Load prob atlas
    prob_nii = nib.load(os.path.join(atlas_dir, 'prob_atlas.nii.gz'))
    # try:
    #     prob_nii = nib.load(os.path.join(atlas_dir), 'prob_atlas.nii.gz')
    # except:
    #     print('* Could not open probabilistic atlas')
    #     return

    prob_atlas = prob_nii.get_data()

    # Loop over labels
    for l in range(0, prob_atlas.shape[3]):

        print('  Creating tryptic for label %03d' % l)

        p = prob_atlas[:,:,:,l]

        # Find minimum isotropic bounding box of p > 0.01 region
        bb = isobb(p > 0.01)

        # Create tryptic of ROI central slices
        tryptic = central_slices(p, isobb(p > 0.01))

        # Write PNG
        tryptic_fname = os.path.join(atlas_dir, 'report', 'tryptic_%03d.png' % l)
        imsave(tryptic_fname, tryptic)


def isobb(mask):
    """
    Determine minimum isotropic bounding box

    Parameters
    ----------
    mask: 3D boolean array

    Returns
    -------
    sx, sy, sz: isotropic bounding box slices
    """

    # Locate objects within mask
    bb = find_objects(mask)

    # Should be only one object, but just in case
    n_obj = len(bb)
    if n_obj > 1:
        print('+ Strange. Found %d objects in label' % n_obj)
        print('+ Using first object')

    # Parse slices for center of mass and maximum length
    sx, sy, sz = bb[0]
    x0 = np.ceil(np.mean([sx.start, sx.stop])).astype(int)
    y0 = np.ceil(np.mean([sy.start, sy.stop])).astype(int)
    z0 = np.ceil(np.mean([sz.start, sz.stop])).astype(int)

    xw = sx.stop - sx.start
    yw = sy.stop - sy.start
    zw = sz.stop - sz.start

    w = np.max([xw, yw, zw]).astype(int)

    return x0, y0, z0, w


def central_slices(p, roi):
    """

    Parameters
    ----------
    img: 3D numpy array
    roi: tuple containing (x0, y0, z0, w) for ROI

    Returns
    -------
    pp: numpy tryptic of central slices
    """

    # Unpack ROI
    x0, y0, z0, w = roi

    # Define central slices
    xx = slice(x0 - w, x0 + w, 1)
    yy = slice(y0 - w, y0 + w, 1)
    zz = slice(z0 - w, z0 + w, 1)

    # Extract slices
    p_xy = p[xx,yy,z0]
    p_xz = p[xx,y0,zz]
    p_yz = p[x0,yy,zz]

    # Create tryptic
    pp = np.concatenate([p_xy, p_xz, p_yz], axis=1)

    return pp


def intra_observer_reports(report_dir, intra_stats):
    """
    Generate intra-observer reports (one per observer)

    Parameters
    ----------
    report_dir: report directory path
    intra_stats: numpy array of intra-observer stats (see load_metrics())

    Returns
    -------

    """

    # Parse stats array
    # Columns: LabelName, LabelNo, Observer, TemplateA, TemplateB, Dice, Hausdorff

    label_names, label_nos, observers, tmp_a, tmp_b, dice, haus = np.transpose(intra_stats)

    # Eliminate duplicates
    label_names = np.unique(label_names)
    label_nos = np.unique(label_nos)
    observers = np.unique(observers)
    tmp_a = np.unique(tmp_a)

    # Count unique labels, observers, templates
    n_lab = len(label_names)
    n_obs = len(observers)
    n_tmp = len(tmp_a)

    # Reshape Dice and Hausdorff matrices
    dice_matrix = dice.reshape(n_lab, n_obs, n_tmp, n_tmp)
    haus_matrix = haus.reshape(n_lab, n_obs, n_tmp, n_tmp)

    # Setup Jinja2 template
    template_loader = jinja2.FileSystemLoader(searchpath="/Users/jmt/GitHub/atlaskit")
    template_env = jinja2.Environment(loader=template_loader)
    template_fname = "atlas_intra_observer.jinja"
    template = template_env.get_template(template_fname)

    # Loop over observers
    for obs in observers:

        # Extract intra similarities for this observer
        this_dice = dice_matrix[:, obs, :, :]
        this_haus = haus_matrix[:, obs, :, :]

        # Intra similarity image filenames
        intra_dice_img = os.path.join(report_dir, "intra_dice_obs_%02d.png" % obs)
        intra_haus_img = os.path.join(report_dir, "intra_haus_obs_%02d.png" % obs)

        # Construct images
        intra_similarity_image(intra_dice_img, this_dice)
        intra_similarity_image(intra_haus_img, this_haus)

        # Template variables
        template_vars = {"observer": obs,
                         "intra_dice_img": intra_dice_img,
                         "intra_haus_img": intra_haus_img}

        # Render page
        output_text = template.render(template_vars)

        # Write report for this observer
        obs_html = os.path.join(report_dir, "intra_observer_%02d.html" % obs)
        with open(obs_html, "w") as f:
            f.write(output_text)


def inter_observer_report(report_dir, inter_stats, label_key):
    """

    Parameters
    ----------
    report_dir
    inter_stats
    label_key

    Returns
    -------

    """

    # Unique labels (column 1, trueidx)
    label_nos = np.int32(np.unique(inter_stats[:,1]))
    n_lab = len(label_nos)

    # Inter stats columns: idx, trueidx, tmp, obsA, obsB, dice, hausdorff
    n_tmp = len(np.unique(inter_stats[:,2]))
    n_obs = len(np.unique(inter_stats[:,3]))

    # Reshape dice and hausdorff matrices
    dice_matrix = inter_stats[:, 5].reshape(n_lab, n_tmp, n_obs, n_obs)
    haus_matrix = inter_stats[:, 6].reshape(n_lab, n_tmp, n_obs, n_obs)

    # Setup template
    template_loader = jinja2.FileSystemLoader(searchpath="/Users/jmt/GitHub/atlaskit")
    template_env = jinja2.Environment(loader=template_loader)
    template_fname = "atlas_intra_observer.jinja"
    template = template_env.get_template(template_fname)

    # Inter similarity image filenames
    inter_dice_img = os.path.join(report_dir, "inter_dice.png")
    inter_haus_img = os.path.join(report_dir, "inter_haus.png")

    # Construct images
    inter_similarity_image(inter_dice_img, dice_matrix, label_nos, label_key)
    inter_similarity_image(inter_haus_img, haus_matrix, label_nos, label_key)

    # Template variables
    template_vars = {"inter_dice_img": inter_dice_img,
                     "inter_haus_img": inter_haus_img}

    # Finally, process the template to produce our final text.
    output_text = template.render(template_vars)

    # Write page to report directory
    obs_html = os.path.join(report_dir, "inter_observer.html")
    with open(obs_html, "w") as f:
        f.write(output_text)


def intra_similarity_image(intra_img, intra_mat, label_nos, label_key):
    """
    Create intra-observer similarity matrix image

    Parameters
    ----------
    intra_img: output image path
    intra_mat: intra-obs similarity coefficient matrix (n_label x n_tmp x n_tmp)

    Returns
    -------

    """

    # Loop over labels
    for l in label_nos:

        label_name = get_label_name(l, label_key)

        plt.matshow(m)
        plt.title(label_name)
        plt.imsave(intra_img)


def inter_similarity_image(inter_img, inter_mat):
    """
    Create inter-observer similarity matrix image

    Parameters
    ----------
    inter_img: output image path
    inter_mat: inter-obs similarity coefficient matrix (n_label x n_obs x n_obs)

    Returns
    -------

    """


def load_metrics(atlas_dir):
    """
    Parse similarity metrics from CSV file

    Parameters
    ----------
    fname : CSV filename to parse

    Returns
    -------
    m : numpy array containing label, observer and template indices and metrics

    """

    # For intra-observer file, columns are idx, true idx, obs, tmpA, tmpB, dice, hausdorff

    intra_csv = os.path.join(atlas_dir, 'intra_observer_metrics.csv')
    with open(intra_csv, "r") as f:
        reader = csv.reader(f)
        l = list(reader)

    # Parse metrics list

    m = np.array(l[1:])

    return intra_stats, inter_stats


def get_template_ids(label_dir, obs):

    obs_dir = os.path.join(label_dir, obs)
    if os.path.isdir(obs_dir):
        ims = glob(os.path_join(obs_dir, '*.nii.gz'))
        print(ims)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()