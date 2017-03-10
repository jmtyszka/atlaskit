#!/usr/bin/env python3
"""
Warp a probabistic atlas onto a brain containing a lesion.
Requires a T1w image of the lesioned brain and a binary lesion mask, both in the same space.
This command outputs lesion volumes within each probabilistic atlas label

Usage
----
atlas2lesion.py -i <lesion T1w> -m <lesion mask> -t <atlas T1w template> -p <probabilistic atlas>
atlas2lesion.py -h

Example
----
>>> prob_label_volumes.py *_probs.nii.gz

Authors
----
Mike Tyszka, Caltech Brain Imaging Center

Dates
----
2016-09-26 JMT From scratch

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
2016 California Institute of Technology.
"""

__version__ = '0.1.0'

import os
import sys
import argparse
import nibabel as nib
import numpy as np
import nipype


def main():

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Atlas-based lesion volumetrics')
    parser.add_argument('prob_files', type=str, nargs='+', help="List of 4D prob label images")

    # Parse command line arguments
    args = parser.parse_args()
    prob_files = args.prob_files

    # Clean exit
    sys.exit(0)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
