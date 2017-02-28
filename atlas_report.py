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


def main():

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Create labeling report for a probabilistic atlas')
    parser.add_argument('-a','--atlasdir', help='Directory containing probabilistic atlas')
    parser.add_argument('-k', '--key', help='ITK-SNAP label key file [labels.txt]')

    # Parse command line arguments
    args = parser.parse_args()

    if args.atlasdir:
        atlas_dir = args.atlasdir
    else:
        atlas_dir = os.path.realpath(os.getcwd())

    print('Atlas directory : %s' % atlas_dir)

    # Check for label directory existence
    if not os.path.isdir(atlas_dir):
        print('Atlas directory does not exist (%s) - exiting' % atlas_dir)
        sys.exit(1)

    # Load and parse label key if provided
    if args.key:
        label_key = load_key(os.path.join(atlas_dir, args.key))
    else:
        label_key = []

    # Create report directory within atlas directory
    report_dir = os.path.join(atlas_dir, 'report')
    if not os.path.isdir(report_dir):
        os.mkdir(report_dir)
    print('Report directory : %s' % report_dir)

    # Inter/intra-observer metrics CSV files
    inter_metrics_csv = os.path.join(atlas_dir, 'inter_observer_metrics.csv')
    intra_metrics_csv = os.path.join(atlas_dir, 'intra_observer_metrics.csv')

    print('Loading similarity metrics')
    intra_stats = load_metrics(intra_metrics_csv)
    inter_stats = load_metrics(inter_metrics_csv)

    # Summary report page
    print('Writing report summary page')
    summary_report(report_dir, intra_stats, inter_stats)

    # Observer report pages
    # observer_reports(intra_stats, inter_stats)

    # Label report pages
    # label_reports()

    # Clean exit
    sys.exit(0)


def summary_report(report_dir, intra_stats, inter_stats):

    template_loader = jinja2.FileSystemLoader(searchpath="/Users/jmt/GitHub/atlaskit")
    template_env = jinja2.Environment(loader=template_loader)
    template_fname = "atlas_summary.jinja"
    template = template_env.get_template(template_fname)

    # Specify any input variables to the template as a dictionary.
    templateVars = {"title": "Test Example",
                    "description": "A simple inquiry of function."}

    # Finally, process the template to produce our final text.
    output_text = template.render(templateVars)

    # Write page to report directory
    with open(os.path.join(report_dir, 'index.html'), "w") as f:
        f.write(output_text)


def load_metrics(fname):
    """
    Parse similarity metrics from CSV file

    Parameters
    ----------
    fname : CSV filename to parse

    Returns
    -------

    """
    with open(fname, "r") as f:
        reader = csv.reader(f)
        l = list(reader)

    m = np.array(l[1:], dtype=np.float)

    dice, haus = m[:,5], m[:,6]

    return dice, haus


def get_template_ids(label_dir, obs):

    obs_dir = os.path.join(label_dir, obs)
    if os.path.isdir(obs_dir):
        ims = glob(os.path_join(obs_dir, '*.nii.gz'))
        print(ims)


def load_key(key_fname):
    '''
    Parse an ITK-SNAP label key file

    Parameters
    ----------
    key_fname

    Returns
    -------

    '''

    # Import key as a data table
    # Note the partially undocumented delim_whitespace flag
    data = pd.read_table(key_fname,
                         comment='#',
                         header=None,
                         names=['Index','R','G','B','A','Vis','Mesh','Name'],
                         delim_whitespace=True)

    return data


def get_label_name(label_idx, label_key):
    '''
    Search label key for label index and return name

    Parameters
    ----------
    label_idx
    label_key

    Returns
    -------

    '''

    label_name = 'Unknown Label'

    for i, idx in enumerate(label_key.Index):
        if label_idx == idx:
            label_name = label_key.Name[i]

    return label_name


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()