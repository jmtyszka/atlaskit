#!/usr/bin/env python3
"""
Create a report of intra and inter-observer atlas label statistics
- requires that atlas.py has been run previously on the labels directory
- generates HTML report pages in subdirectory of atlas directory

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
import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
from scipy.misc import imresize
from scipy.ndimage import find_objects
from datetime import datetime


__version__ = '1.0.0'


def main():

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Create labeling report for a probabilistic atlas')
    parser.add_argument('-a','--atlasdir', required=True, help='Directory containing probabilistic atlas')

    # Parse command line arguments
    args = parser.parse_args()
    atlas_dir = args.atlasdir


    print('')
    print('-----------------------------')
    print('Atlas label similarity report')
    print('-----------------------------')

    # Check for atlas directory existence
    if not os.path.isdir(atlas_dir):
        print('Atlas directory does not exist (%s) - exiting' % atlas_dir)
        sys.exit(1)

    # Create report directory within atlas directory
    report_dir = os.path.join(atlas_dir, 'report')
    if not os.path.isdir(report_dir):
        os.mkdir(report_dir)

    print('Atlas directory  : %s' % atlas_dir)
    print('Report directory : %s' % report_dir)
    print('')

    print('Loading similarity metrics')
    intra_stats, inter_stats = load_metrics(atlas_dir)

    # Intra-observer reports (one per observer)
    print('')
    print('Generating intra-observer report')
    intra_observer_report(report_dir, intra_stats)

    # Inter-observer report
    print('')
    print('Generating inter-observer report')
    inter_observer_report(report_dir, inter_stats)

    # Summary report page
    print('')
    print('Writing report summary page')
    summary_report(atlas_dir, report_dir, intra_stats, inter_stats)

    # Clean exit
    sys.exit(0)


def summary_report(atlas_dir, report_dir, intra_metrics, inter_metrics):
    """
    Summary report for the entire atlas
    - maximum probability projections for all labels

    Parameters
    ----------
    atlas_dir: atlas directory path
    report_dir: report directory path
    intra_stats: intra-observer metrics tuple
    inter_stats: inter-observer metrics tuple

    Returns
    -------

    """

    # Setup Jinja2 template
    html_loader = jinja2.FileSystemLoader(searchpath="./")
    html_env = jinja2.Environment(loader=html_loader)
    html_fname = "atlas_summary.jinja"
    html = html_env.get_template(html_fname)

    # Parse metrics tuples
    label_names, label_nos, observers, templates, intra_dice, intra_haus = intra_metrics
    _, _, _, _, inter_dice, inter_haus = inter_metrics

    # Determine montage size from number of labels
    ncols = 4
    nrows = np.ceil(len(label_names)/ncols).astype(int)

    # Create tryptic overlays through CoM of label
    print('  Generating maximum probability projections for each label')
    mpp_fname = maxprob_projections(atlas_dir, report_dir, label_names, nrows, ncols)

    # Init observer stats list to pass to HTML template
    stats = []

    print('  Generating similarity summary tables')

    # Construct list of results for each label
    for obs in observers:

        print('    Observer %d' % obs)

        obs_stats = []

        for ll, label_name in enumerate(label_names):

            this_intra_dice = intra_dice[ll, obs, :, :]
            this_intra_haus = intra_haus[ll, obs, :, :]

            intra_dice_mean = mean_str(this_intra_dice)
            intra_haus_mean = mean_str(this_intra_haus)

            # Find unfinished template labels
            # Search for NaNs on leading diagonals in intra dice data
            unfinished = str(np.where(np.isnan(np.diagonal(this_intra_dice)))[0])

            label_dict = dict([("label_name", label_name),
                               ("label_no", label_nos[ll]),
                               ("intra_dice_mean", intra_dice_mean),
                               ("intra_haus_mean", intra_haus_mean),
                               ("unfinished", unfinished)])

            obs_stats.append(label_dict)

        stats.append(obs_stats)

    # Template variables
    template_vars = {"mpp_fname": mpp_fname,
                     "stats": stats,
                     "report_time": datetime.now().strftime('%Y-%m-%d %H:%M')}

    # Finally, process the template to produce our final text.
    output_text = html.render(template_vars)

    # Write page to report directory
    with open(os.path.join(os.path.join(atlas_dir), 'report', 'index.html'), "w") as f:
        f.write(output_text)


def intra_observer_report(report_dir, intra_metrics):
    """
    Generate intra-observer report for each observer

    Parameters
    ----------
    report_dir: report directory path
    intra_metrics: tuple containing labelNames, labelNos, observers, templates, dice and haussdorff metrics

    Returns
    -------

    """

    # Setup Jinja2 template
    html_loader = jinja2.FileSystemLoader(searchpath="./")
    html_env = jinja2.Environment(loader=html_loader)
    html_fname = "atlas_intra_observer.jinja"
    html = html_env.get_template(html_fname)

    # Parse metrics tuple
    label_names, label_nos, observers, templates, dice, haus = intra_metrics

    # Determine montage size from number of labels
    ncols = 4
    nrows = np.ceil(len(label_names)/ncols).astype(int)

    # Metric limits
    dlims = 0.0, 1.0
    hlims = 0.0, 10.0

    # Create similarity figures over all labels and observers
    intra_dice_imgs = similarity_figure(dice, observers, "Observer %d Dice Coefficient", "intra_obs_%0d_dice.png",
                                        report_dir, label_names, dlims, nrows, ncols, 0.0)
    intra_haus_imgs = similarity_figure(haus, observers, "Observer %d Hausdorff Distance (mm)", "intra_obs_%0d_haus.png",
                                        report_dir, label_names, hlims, nrows, ncols, 1e6)

    # Composite all images into a single dictionary list
    intra_imgs = []
    for i, dimg in enumerate(intra_dice_imgs):
        himg = intra_haus_imgs[i]
        intra_imgs.append(dict(dimg=dimg, himg=himg))

    # Template variables
    html_vars = {"intra_imgs": intra_imgs,
                 "report_time": datetime.now().strftime('%Y-%m-%d %H:%M')}

    # Render page
    html_text = html.render(html_vars)

    # Write report
    obs_html = os.path.join(report_dir, "intra_report.html")
    with open(obs_html, "w") as f:
        f.write(html_text)


def inter_observer_report(report_dir, inter_metrics):
    """
    Generate inter-observer report for each template

    Parameters
    ----------
    report_dir: report directory path
    inter_metrics: tuple containing labelNames, labelNos, observers, templates, dice and haussdorff metrics

    Returns
    -------

    """

    # Setup Jinja2 template
    html_loader = jinja2.FileSystemLoader(searchpath="./")
    html_env = jinja2.Environment(loader=html_loader)
    html_fname = "atlas_inter_observer.jinja"
    html = html_env.get_template(html_fname)

    # Parse metrics tuple
    label_names, label_nos, observers, templates, dice, haus = inter_metrics

    # Determine subplot matrix dimensions from number of labels
    # nrows x ncols where nrows = ceil(n_labels/ncols)
    ncols = 4
    nrows = np.ceil(len(label_names)/ncols).astype(int)

    # Metric limits
    dlims = 0.0, 1.0
    hlims = 0.0, 10.0

    # Create similarity figures over all labels and observers
    inter_dice_imgs = similarity_figure(dice, templates, "Template %d : Dice Coefficient", "inter_tmp_%0d_dice.png",
                                        report_dir, label_names, dlims, nrows, ncols, 0.0)
    inter_haus_imgs = similarity_figure(haus, templates, "Template %d Hausdorff Distance (mm)", "inter_tmp_%0d_haus.png",
                                        report_dir, label_names, hlims, nrows, ncols, 1e6)

    # Composite all images into a single dictionary list
    inter_imgs = []
    for i, dimg in enumerate(inter_dice_imgs):
        himg = inter_haus_imgs[i]
        inter_imgs.append(dict(dimg=dimg, himg=himg))

    # Template variables
    html_vars = {"inter_imgs": inter_imgs,
                 "report_time": datetime.now().strftime('%Y-%m-%d %H:%M')}

    # Render page
    html_text = html.render(html_vars)

    # Write report
    obs_html = os.path.join(report_dir, "inter_report.html")
    with open(obs_html, "w") as f:
        f.write(html_text)


def maxprob_projections(atlas_dir, report_dir, label_names, nrows, ncols):
    """
    Construct an array of maximum probablity projections through each label
    over all observers and templates

    Parameters
    ----------
    atlas_dir: atlas directory path
    report_dir: report directory path
    label_names: list of unique label names (in label number order)
    nrows, ncols: figure matrix size

    Returns
    -------
    mpp_png: maxprob projection PNG filename
    """

    # Probability threshold for minimum BB
    p_thresh = 0.25

    # Load prob atlas
    prob_nii = nib.load(os.path.join(atlas_dir, 'prob_atlas.nii.gz'))
    prob_atlas = prob_nii.get_data()

    # Create figure with subplot array
    fig, axs = plt.subplots(nrows, ncols, figsize=(8,4))
    axs = np.array(axs).reshape(-1)

    # Loop over axes
    for aa, ax in enumerate(axs):

        if aa < len(label_names):

            print('    %s' % label_names[aa])

            # Current prob label
            p = prob_atlas[:, :, :, aa]

            # Create tryptic of central slices through ROI defined by p > p_thresh
            tryptic = central_slices(p, isobb(p > p_thresh))

            ax.pcolor(tryptic)
            ax.set_title(label_names[aa], fontsize=8)

        else:
            ax.axis('off')

        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.set(adjustable='box-forced', aspect='equal')
        ax.set(aspect='equal')

    # Tidy up spacing
    plt.tight_layout()

    # Save figure to PNG
    mpp_fname = 'mpp.png'
    plt.savefig(os.path.join(report_dir, mpp_fname))

    # Clean up
    plt.close(fig)

    return mpp_fname


def similarity_figure(metric, inds, tfmt, ffmt, report_dir, label_names, mlims, nrows, ncols, nansub=0.0):
    """
    Plot an array of similarity matrix figures

    Parameters
    ----------
    metric: similarity metric array to plot
    inds: index iterator (observers or templates)
    tfmt: title format for each figure
    ffmt: file format string
    report_dir: report directory
    label_names: list of label names
    mlims: scale limits for metric
    nrows: plot grid rows
    ncols: plot grid columns
    nansub: value to replace NaNs in data

    Returns
    -------
    img_list: list of image filenames
    """

    # Init image file list
    img_list = []

    # Loop over indices (observers or templates)
    for ii in inds:

        # Create figure with subplot array
        fig, axs = plt.subplots(nrows, ncols)
        axs = np.array(axs).reshape(-1)

        # Construct subplot matrix
        mm = metric[:, ii, :, :]
        for aa, ax in enumerate(axs):

            if aa < len(label_names):
                mmaa = np.flipud(mm[aa, :, :]).copy()
                mmaa[np.isnan(mmaa)] = nansub
                im = ax.pcolor(mmaa, vmin=mlims[0], vmax=mlims[1])
                ax.set_title(label_names[aa], fontsize=8)
            else:
                ax.axis('off')

            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            ax.set(adjustable='box-forced', aspect='equal')

        # Tidy up spacing
        plt.tight_layout()

        # Make space for title and colorbar
        fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8)
        plt.suptitle(tfmt % ii, x=0.5, y=0.99)
        cax = fig.add_axes([0.85, 0.1, 0.05, 0.8]) # [x0, y0, w, h]
        fig.colorbar(im, cax=cax)

        # Save figure to PNG
        fname = ffmt % ii
        print('  %s' % fname)
        img_list.append(fname)
        plt.savefig(os.path.join(report_dir, fname))

        # Clean up
        plt.close(fig)

    return img_list


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

    if any(bb):

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

    else:

        x0, y0, z0, w = 0, 0, 0, -1

    return x0, y0, z0, w


def central_slices(p, roi):
    """
    Create tryptic of central slices from ROI

    Parameters
    ----------
    img: 3D numpy array
    roi: tuple containing (x0, y0, z0, w) for ROI

    Returns
    -------
    pp: numpy tryptic of central slices
    """

    # Base output image dimensions
    base_dims = 64, 3*64

    # Unpack ROI
    x0, y0, z0, w = roi

    if w > 0:

        # Define central slices
        xx = slice(x0 - w, x0 + w, 1)
        yy = slice(y0 - w, y0 + w, 1)
        zz = slice(z0 - w, z0 + w, 1)

        # Extract slices
        p_xy = p[xx,yy,z0]
        p_xz = p[xx,y0,zz]
        p_yz = p[x0,yy,zz]

        # Create horizontal tryptic
        pp = np.hstack([p_xy, p_xz, p_yz])

        # Resize to base size
        pp = imresize(pp, base_dims, interp='bicubic')

    else:

        # Blank tryptic
        pp = np.zeros(base_dims)

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

    # Find unique label numbers with initial row indices for each
    label_nos, idx = np.unique(m['labelNo'], return_index=True)

    # Extract corresponding label names to unique label numbers
    label_names = m['labelName'][idx].astype(str)

    # Unique template and observer lists can be sorted as usual
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

    # Find unique label numbers with initial row indices for each
    label_nos, idx = np.unique(m['labelNo'], return_index=True)

    # Extract corresponding label names to unique label numbers
    label_names = m['labelName'][idx].astype(str)

    # Unique template and observer lists can be sorted as usual
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


def mean_str(x):
    """
    Create formated string containing mean of x[]
    Handle NaNs

    Parameters
    ----------
    x: numpy float array

    Returns
    -------
    xms: formatted mean string

    """

    n_nans = np.sum(np.isnan(x))

    if n_nans == x.size:
        xms = "-"
    else:
        xms = "%0.3f" % np.nanmean(x)

    return xms

# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
