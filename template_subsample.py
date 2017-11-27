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
import re
import argparse
import random as rnd
import numpy as np
from nipype.interfaces.ants import AverageImages


def main():

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Generate participant subsamples for midspace template variants')
    parser.add_argument('-b', '--bidsdir', default='.',
                        help='BIDS study directory containing source and derivatives folders [.]')
    parser.add_argument('-n','--nvariants', default=1,
                        help='Number of template variants [1]')

    # Parse command line arguments
    args = parser.parse_args()

    bids_dir = args.bidsdir
    nv = args.nvariants

    # Derivatives directory
    der_dir = os.path.join(bids_dir, 'derivatives')

    # Midspace directory (containing midspace transformed T1w images)
    midspace_dir = os.path.join(der_dir, 'midspace')

    # Output template directory
    template_dir = os.path.join(der_dir, 'templates')
    if not os.path.isdir(template_dir):
        os.makedirs(template_dir, exist_ok=True)

    print('BIDS Study Directory : %s' % bids_dir)
    print('Midspace Directory   : %s' % midspace_dir)
    print('Templates Directory  : %s' % template_dir)

    # Construct list of all individual images in midspace directory
    img_list = [f for f in os.listdir(midspace_dir) if re.search(r'[0-9]Warp\.nii\.gz$', f)]

    # Total number of images
    n_img = len(img_list)

    if n_img < 1:
        print('* No individual midspace images found in derivatives/midspace directory')
        print('* Exiting')
        sys.exit(1)

    # Check for odd number of images
    is_odd = bool(n_img % 2)

    hn = int(n_img / 2)

    # Standard random seed
    rnd.seed(1966)

    # Loop over requested number of variants
    # Generates 2n participant lists
    for vc in np.arange(nv):

        print('Variant %d' % (vc+1))

        # Generate randomized image index list
        idx = np.arange(0, n_img)
        rnd.shuffle(idx)

        # Discard first index from randomized list if odd number of images
        if is_odd:
            print('* Odd number of images - dropping %s' % img_list[0])
            idx = np.delete(idx, 0)

        # Divide list into A and B sets
        idx_A = idx[0:hn]
        idx_B = idx[hn:]

        tname_A = os.path.join(template_dir, 'template_%03d_A.nii.gz' % vc)
        tname_B = os.path.join(template_dir, 'template_%03d_A.nii.gz' % vc)

        # Average images from sample A
        average_images(img_list, idx_A, tname_A)

        # Average images from conjugate sample B
        average_images(img_list, idx_B, tname_B)


def average_images(img_list, idx, tname):
    """

    Parameters
    ----------
    img_list: list of strings
    idx: numpy array
    tname: string

    Returns
    -------

    """

    # Init subsampled image list
    sample_imgs = []

    # Save image list to sidecar .txt file
    tname_txt = tname.replace('nii.gz','txt')
    with open(tname_txt, 'w') as fd:
        for ii in idx:
            sample_imgs.append(img_list[ii])
            fd.write(img_list[ii] + '\n')

    # Call ANTs AverageImages command using nipype
    avg = AverageImages()
    avg.inputs.dimension = 3
    avg.inputs.output_average_image = tname
    avg.inputs.normalize = False
    avg.inputs.images = sample_imgs
    avg.run()

    return img_list


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()

