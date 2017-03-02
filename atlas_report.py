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
from glob import glob
import nibabel as nib
import multiprocessing as mp
import matplotlib.pyplot as plt


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

    # Load the label key from the atlas directory as a data frame
    label_key = load_key(os.path.join(atlas_dir, 'labels.txt'))

    # Create report directory within atlas directory
    report_dir = os.path.join(atlas_dir, 'report')
    if not os.path.isdir(report_dir):
        os.mkdir(report_dir)
    print('Report directory : %s' % report_dir)

    print('Loading similarity metrics')
    inter_metrics_csv = os.path.join(atlas_dir, 'inter_observer_metrics.csv')
    intra_metrics_csv = os.path.join(atlas_dir, 'intra_observer_metrics.csv')
    intra_stats = load_metrics(intra_metrics_csv)
    inter_stats = load_metrics(inter_metrics_csv)

    # Determine label numbers present (column 1, trueidx)
    label_nos = np.int32(np.unique(intra_stats[:,1]))

    # Intra-observer reports (one per observer)
    print('Generating intra-observer report')
    intra_observer_reports(report_dir, intra_stats, label_nos, label_key)

    # Inter-observer report
    print('Generating inter-observer report')
    inter_observer_report(report_dir, inter_stats, label_nos, label_key)

    # Summary report page
    # print('Writing report summary page')
    # summary_report(report_dir, )

    # Clean exit
    sys.exit(0)


def summary_report(report_dir, intra_stats, inter_stats, label_nos, label_key):
    """

    Parameters
    ----------
    report_dir: report directory path
    intra_stats: numpy array of intra-observer stats (see load_metrics())
    inter_stats: numpy array of inter-observer stats (see load_metrics())

    Returns
    -------

    """

    # Intra stats columns: idx, trueidx, obs, tmpA, tmpB, dice, hausdorff
    n_obs = len(np.unique(intra_stats[:,2]))
    n_tmp = len(np.unique(intra_stats[:,3]))

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
    with open(os.path.join(report_dir, 'index.html'), "w") as f:
        f.write(output_text)


def intra_observer_reports(report_dir, intra_stats, label_nos, label_key):
    """
    Generate intra-observer reports (one per observer)

    Parameters
    ----------
    report_dir: report directory path
    intra_stats: numpy array of intra-observer stats (see load_metrics())
    label_nos: list of unique labels in the atlas
    label_key: label key dictionary

    Returns
    -------

    """

    # Intra stats columns: idx, trueidx, obs, tmpA, tmpB, dice, hausdorff
    n_lab = len(label_nos)
    n_obs = len(np.unique(intra_stats[:,2]))
    n_tmp = len(np.unique(intra_stats[:,3]))

    # Reshape Dice and Hausdorff matrices
    dice_matrix = intra_stats[:,5].reshape(n_lab, n_obs, n_tmp, n_tmp)
    haus_matrix = intra_stats[:,6].reshape(n_lab, n_obs, n_tmp, n_tmp)

    # Setup template
    template_loader = jinja2.FileSystemLoader(searchpath="/Users/jmt/GitHub/atlaskit")
    template_env = jinja2.Environment(loader=template_loader)
    template_fname = "atlas_intra_observer.jinja"
    template = template_env.get_template(template_fname)

    # Loop over observers
    for obs in range(0, n_obs):

        # Extract intra similarities for this observer
        this_dice = dice_matrix[:,obs,:,:]
        this_haus = haus_matrix[:,obs,:,:]

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
    unique_labels = np.int32(np.unique(inter_stats[:,1]))
    n_lab = len(unique_labels)

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
    inter_similarity_image(inter_haus_img, haus_matrix)

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


def load_metrics(fname):
    """
    Parse similarity metrics from CSV file

    Parameters
    ----------
    fname : CSV filename to parse

    Returns
    -------
    m : numpy array containing label, observer and template indices and metrics
    For intra-observer file, columns are idx, true idx, obs, tmpA, tmpB, dice, hausdorff
    For inter-observer file, columns are idx, true idx, tmp, obsA, obsB, dice, hausdorff

    """
    with open(fname, "r") as f:
        reader = csv.reader(f)
        l = list(reader)

    m = np.array(l[1:], dtype=np.float)

    return m


def get_template_ids(label_dir, obs):

    obs_dir = os.path.join(label_dir, obs)
    if os.path.isdir(obs_dir):
        ims = glob(os.path_join(obs_dir, '*.nii.gz'))
        print(ims)


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

    # Import key as a data table
    # Note the partially undocumented delim_whitespace flag
    key = pd.read_table(key_fname,
                         comment='#',
                         header=None,
                         names=['Index','R','G','B','A','Vis','Mesh','Name'],
                         delim_whitespace=True)

    return key


def get_label_name(label_idx, label_key):
    """
    Search label key for label index and return name

    Parameters
    ----------
    label_idx
    label_key

    Returns
    -------

    """

    label_name = 'Unknown Label'

    for i, idx in enumerate(label_key.Index):
        if label_idx == idx:
            label_name = label_key.Name[i]

    return label_name


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()