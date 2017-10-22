#!/usr/bin/env python3
"""
Output volumes of each probabilistic label in an atlas in microliters.
Accepts multiple 4D prob atlas files

Usage
----
prob_label_volumes.py <atlas_file>
prob_label_volumes.py -h

Example
----
>>> prob_label_volumes.py *_probs.nii.gz

Authors
----
Mike Tyszka, Caltech Brain Imaging Center

Dates
----
2015-05-31 JMT Adapt from label_volumes.py

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
2015 California Institute of Technology.
"""

__version__ = '0.1.0'

import os
import sys
import argparse
import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt


def load_labels(filename):
    ''' Load labels from file, which was used as label key in inkscape '''

    keys = []
    
    label_file = open(filename, 'r')
    entries = label_file.readlines()
    label_file.close()

    for e in range(0, len(entries)):
        entry = entries[e]
        start = entry.find('"') + 1
        end = entry.find('"', start)
        keys.append(entry[start:end])
    return keys


def create_histogram(p, keys, nrows=4, ncols=4, fontsize=14, img_fname='/tmp/prob_atlas_hist.png'):
    ''' for each label in probabilistic atlas, create a histogram of probabilities '''
    
    fig, axs = plt.subplots(nrows, ncols)
    axs = np.array(axs).reshape(-1)

    im = []

    for aa, ax in enumerate(axs):

        if aa < len(keys):
            v = p[:,:,:,aa].ravel()
            v = v[v>0]
            im = ax.hist(v, bins=50, cumulative=True, density=True, range=(0,1))
            ax.set_title(keys[aa], fontsize=fontsize)
        else:
            ax.axis('off')

        if aa > 0:
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
        else:
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(fontsize)        
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(fontsize)        
            
    # Tidy up spacing
    plt.tight_layout()

    print('  Saving image to %s' % img_fname)
    plt.savefig(os.path.join(img_fname), bbox_inches='tight')

def print_vol_info(p_file, keys, latex=False):

    if latex:
        print("\\begin{tabular}{l l}\n\hline\\\\\nLabel&Vol (muL)\\\\\n\hline\\\\")
    # Load the source atlas image
    p_nii = nib.load(p_file)
    p = p_nii.get_data()
    nd = p.ndim
        
    # Atlas voxel volume in mm^3 (microliters)
    atlas_vox_vol_ul = np.array(p_nii.header.get_zooms()).prod()
                   
    nx,ny,nz,nt = p.shape

    for t in range(0,nt):
        label = keys[t]
        V = np.sum(p[:,:,:,t]) * atlas_vox_vol_ul
        if latex:
            print('%s & %0.3f\\\\' % (label.replace('_','\_'), V))
        else:
            print('%s:\t%0.3f' % (label, V))

    if latex:
        print("\hline\\\\\\end{tabular}\n")
    # Final newline
    print
    
    return p
    
def main():
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Probabilistic label volumes in microliters')
    parser.add_argument('-d', '--atlas_dir', required=False, help="Director of probabilistic atlas and label file")
    parser.add_argument('-f', '--p_file', required=False, nargs='+', help="List of 4D prob label images")
    parser.add_argument('-l', '--l_file', required=False, nargs='+', help='List of (itksnap) label files')
    parser.add_argument('--latex', dest='latex', action='store_true')

    # Parse command line arguments
    args = parser.parse_args()
    atlas_dir = args.atlas_dir
    if not atlas_dir is None:
        print("atlas_dir: %s" % atlas_dir)
        p_file = os.path.join(atlas_dir, 'prob_atlas.nii.gz')
        l_file = os.path.join(atlas_dir, 'labels.txt')
    else:
        print("atlas_dir: not specified")
        p_file = args.p_file
        l_file = args.l_file

    if os.path.exists(p_file):
        print("Probabilistic atlas file: %s" % p_file)
    if os.path.exists(l_file):
        print("Label file: %s" % l_file)        

    latex = args.latex
    
    keys = load_labels(l_file)
        
    # Force absolute path
    p_file = os.path.abspath(p_file)
                
    p = print_vol_info(p_file, keys, latex=latex)
    
    # create histogram of label probabilities for each label
    create_histogram(p, keys)
    
    # Clean exit
    sys.exit(0)

    
# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
