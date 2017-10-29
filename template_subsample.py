#!/usr/bin/env python
"""
Create a set of subsamples from a BIDS participant list
Subsamples are used to generate template variants for a midspace

Authors
----
Mike Tyszka, Caltech Brain Imaging Center

Dates
----
2017-10-25 JMT From scratch

MIT License

Copyright (c) 2017 Mike Tyszka

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

__version__ = '0.1.0'

import os
import sys
import argparse
import pandas as pd


def main():

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Generate participant subsamples for midspace template variants')
    parser.add_argument('-i', '--indir', default='source', help='BIDS source directory [source]')
    parser.add_argument('-o', '--outdir', default='derivatives', help='BIDS derivatives directory [derivatives]')
    parser.add_argument('-n','--nvariants', default=1, help='Number of template variants [1]')

    # Parse command line arguments
    args = parser.parse_args()

    src_dir = args.indir
    der_dir = args.outdir
    n = args.nvariants

    # Load participant list
    pids_tsv = os.path.join(src_dir, 'participants.tsv')
    pids = pd.read_table(pids_tsv)


    # If number of participants is odd, drop last participant from list

    # Generate template variant pairs

    # Write SID list to derivatives folder


def bids_init(bids_root_dir):
    """
    Initialize root BIDS directory
    :param bids_root_dir: root BIDS directory
    :return participants_fd: participant TSV file descriptor
    """

    # Create template participant TSV file in BIDS root directory
    parts_tsv = os.path.join(bids_root_dir, 'participants.tsv')
    participants_fd = open(parts_tsv, 'w')
    participants_fd.write('participant_id\tsex\tage\n')

    # Create template JSON dataset description
    datadesc_json = os.path.join(bids_root_dir, 'dataset_description.json')
    meta_dict = dict({'BIDSVersion': "1.0.0",
                      'License': "This data is made available under the Creative Commons BY-SA 4.0 International License.",
                      'Name': "The dataset name goes here",
                      'ReferencesAndLinks': "References and links for this dataset go here"})

    # Write JSON file
    bids_write_json(datadesc_json, meta_dict)

    return participants_fd


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
