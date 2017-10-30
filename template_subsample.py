#!/usr/bin/env python
"""
Create randomly subsampled midspace template variants
- Requires BIDS source directory with participants.tsv
- Subsamples are used to generate template variants lists for a midspace
- If odd number of participants, drop final participant
- Sample half of participants randomly for one template
- The remaining half creates a complementary template
- Outputs two lists of participant IDs per requested variant

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
import random as rnd
import numpy as np


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
    nv = args.nvariants

    # Load participant list
    pids_tsv = os.path.join(src_dir, 'participants.tsv')
    print('Reading participant list from %s' % pids_tsv)
    pids = pd.read_table(pids_tsv)

    # If number of participants is odd, drop last participant from list
    n = pids.shape[0]
    if bool(n & 1):
        print('* Odd number of participants (%d)' % n)
        pids = pids[:-1]
        n = pids.shape[0]
        print('* Retaining first %d participants' % n)

    # Standard random seed
    rnd.seed(1966)

    # Loop over requested number of variants
    # Generates 2n participant lists
    for vc in np.arange(nv):

        print('Variant %d' % (vc+1))

        # Generate template variant pairs
        idx = np.arange(0, n)
        rnd.shuffle(idx)

        # Divide into A and B complementary samples
        hn = int(n/2)
        idx_A = idx[0:hn]
        idx_B = idx[hn:]

        # Write PID lists for sample A and B to derivatives folder
        pid_A = pids[idx_A][0]
        pid_B = pids[idx_B][0]

        print(pid_A)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
