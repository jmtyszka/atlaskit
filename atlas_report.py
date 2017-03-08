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
    intra_observer_report(report_dir, intra_stats)

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


def intra_observer_report(report_dir, intra_metrics):
    """
    Generate intra-observer report

    Parameters
    ----------
    report_dir: report directory path
    intra_metrics: tuple containing labelNames, labelNos, observers, templates, dice and haussdorff metrics

    Returns
    -------

    """

    # Setup Jinja2 template
    html_loader = jinja2.FileSystemLoader(searchpath="/Users/jmt/GitHub/atlaskit")
    html_env = jinja2.Environment(loader=html_loader)
    html_fname = "atlas_intra_observer.jinja"
    html = html_env.get_template(html_fname)

    # Parse metrics tuple
    label_names, label_nos, observers, templates, dice, haus = intra_metrics

    # Init image filename lists
    intra_dice_imgs = []
    intra_haus_imgs = []

    # Loop over observers
    for obs in observers:

        # Extract intra similarities for this observer
        dd = dice[:, obs, :, :]
        hh = haus[:, obs, :, :]

        dice_fstub = os.path.join(report_dir, "intra_dice_obs%02d" % obs)

        for idx, label_name in enumerate(label_names):
            fname = dice_fstub + '_%s.png' % label_name
            plt.imsave(fname, dd[idx, :, :])
            intra_dice_imgs.append(dict(observer=obs, label=label_name, href=fname))

        haus_fstub = os.path.join(report_dir, "intra_haus_obs%02d" % obs)

        for idx, label_name in enumerate(label_names):
            fname = haus_fstub + '_%s.png' % label_name
            plt.imsave(fname, hh[idx, :, :])
            intra_haus_imgs.append(dict(observer=obs, label=label_name, href=fname))

    # Template variables
    html_vars = {"observers": observers,
                 "intra_dice_imgs": intra_dice_imgs,
                 "intra_haus_imgs": intra_haus_imgs}

    # Render page
    html_text = html.render(html_vars)

    # Write report
    obs_html = os.path.join(report_dir, "intra_report.html")
    with open(obs_html, "w") as f:
        f.write(html_text)


def inter_observer_report(report_dir, inter_metrics):
    """
    Generate inter-observer report

    Parameters
    ----------
    report_dir: report directory path
    inter_metrics: tuple containing labelNames, labelNos, observers, templates, dice and haussdorff metrics

    Returns
    -------

    """

    # Setup Jinja2 template
    html_loader = jinja2.FileSystemLoader(searchpath="/Users/jmt/GitHub/atlaskit")
    html_env = jinja2.Environment(loader=html_loader)
    html_fname = "atlas_inter_observer.jinja"
    html = html_env.get_template(html_fname)

    # Parse metrics tuple
    label_names, label_nos, observers, templates, dice, haus = inter_metrics

    # Init image filename lists
    inter_dice_imgs = []
    inter_haus_imgs = []

    # Loop over templates
    for tmp in templates:

        # Extract inter-observer similarities for this template
        dd = dice[:, tmp, :, :]
        hh = haus[:, tmp, :, :]

        dice_fstub = os.path.join(report_dir, "inter_dice_tmp%02d" % tmp)

        for idx, label_name in enumerate(label_names):
            fname = dice_fstub + '_%s.png' % label_name
            plt.imsave(fname, dd[idx, :, :])
            inter_dice_imgs.append(dict(template=tmp, label=label_name, href=fname))

        haus_fstub = os.path.join(report_dir, "inter_haus_tmp%02d" % tmp)

        for idx, label_name in enumerate(label_names):
            fname = haus_fstub + '_%s.png' % label_name
            plt.imsave(fname, hh[idx, :, :])
            inter_haus_imgs.append(dict(template=tmp, label=label_name, href=fname))

    # Template variables
    html_vars = {"observers": observers,
                 "intra_dice_imgs": inter_dice_imgs,
                 "intra_haus_imgs": inter_haus_imgs}

    # Render page
    html_text = html.render(html_vars)

    # Write report
    obs_html = os.path.join(report_dir, "inter_report.html")
    with open(obs_html, "w") as f:
        f.write(html_text)


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

    # ----------------------------
    # Load intra-observer metrics
    # Ignore number of voxels in each label (nA, nB) for now
    intra_csv = os.path.join(atlas_dir, 'intra_observer_metrics.csv')
    m = np.genfromtxt(intra_csv,
                      dtype=None,
                      names=['labelName', 'labelNo', 'observer', 'tmpA', 'tmpB', 'dice', 'haus', 'nA', 'nB'],
                      delimiter=',', skip_header=1)

    # Parse array
    label_names = np.unique(m['labelName']).astype(str) # Convert from byte to unicode
    label_nos = np.unique(m['labelNo'])
    observers = np.unique(m['observer'])
    templates = np.unique(m['tmpA'])

    # Count labels, templates and observers
    nLabels, nTemplates, nObservers = len(label_names), len(templates), len(observers)

    # Cast to float and reshape metrics
    dice = m['dice'].reshape(nLabels, nObservers, nTemplates, nTemplates)
    haus = m['haus'].reshape(nLabels, nObservers, nTemplates, nTemplates)

    # Composite into intra_metrics tuple
    intra_metrics = label_names, label_nos, observers, templates, dice, haus

    #----------------------------
    # Load inter-observer metrics
    # Ignore number of voxels in each label (nA, nB) for now
    inter_csv = os.path.join(atlas_dir, 'inter_observer_metrics.csv')
    m = np.genfromtxt(inter_csv,
                      dtype=[('labelName', 'a32'), ('labelNo', 'u8'),
                             ('template', 'u8'), ('obsA', 'u8'), ('obsB', 'u8'),
                             ('dice', 'f8'), ('haus', 'f8'),
                             ('nA', 'u8'), ('nB', 'u8')],
                      delimiter=',', skip_header=1)

    # Parse array
    label_names = np.unique(m['labelName'])
    label_nos = np.unique(m['labelNo'])
    templates = np.unique(m['template'])
    observers = np.unique(m['obsA'])

    # Count labels, templates and observers
    nLabels, nTemplates, nObservers = len(label_names), len(templates), len(observers)

    # Cast to float and reshape metrics
    dice = m['dice'].reshape(nLabels, nTemplates, nObservers, nObservers)
    haus = m['haus'].reshape(nLabels, nTemplates, nObservers, nObservers)

    # Composite into inter_metrics tuple
    inter_metrics = label_names, label_nos, observers, templates, dice, haus

    return intra_metrics, inter_metrics


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()