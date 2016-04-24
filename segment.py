#!/usr/bin/env python3
"""
Simple intensity clustering segmentation of 3D grayscale image

Usage
----
segment.py -i <Grayscale Nifti image> -o <Segmentation image> [-n segments] [-m kmeans]

Authors
----
Mike Tyszka, Caltech Brain Imaging Center

Dates
----
2016-04-24 JMT From scratch

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

import sys
import argparse
import numpy as np
import nibabel as nib
from sklearn.cluster import KMeans
from scipy.ndimage.filters import median_filter


def main():

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Simple image segmenter')
    parser.add_argument('-i', '--input', required=True, help='Input grayscale image')
    parser.add_argument('-o', '--output', required=True, help='Output segmentation labels')
    parser.add_argument('-m', '--method', required=False, help='Clustering method [KMeans]')
    parser.add_argument('-n', '--nclusters', required=False, help='Number of clusters [3]')

    # Parse command line arguments
    args = parser.parse_args()

    in_file = args.input
    out_file = args.output

    if args.method:
        seg_method = args.method
    else:
        seg_method = 'KMeans'

    if args.nclusters:
        n = args.nclusters
    else:
        n = 3

    # Load grayscale image
    print('Loading grayscale image from %s' % in_file)
    in_nii = nib.load(in_file)
    in_img = in_nii.get_data()

    # Save dimensions
    orig_shape = in_img.shape

    # Flatten image for segmentation - see sklearn documentation
    in_img = in_img.reshape(-1,1,order='F')

    # Intesnity cluster segmentation
    print('Segmenting using %s' % seg_method)
    if seg_method == 'KMeans':
        k_means = KMeans(init='k-means++', n_clusters=n)
        k_means.fit(in_img)
        seg_img = k_means.labels_
    else:
        print('Unknown clustering method : %s' % seg_method)
        sys.exit(1)

    # Restore original image shape
    seg_img = seg_img.reshape(orig_shape, order='F')

    # Isolated voxel removal
    seg_img = median_filter(seg_img, size=(3,3,3))

    # Write segmentation labels
    print('Saving segmentation to %s' % out_file)
    out_nii = nib.Nifti1Image(seg_img, in_nii.get_affine())
    out_nii.to_filename(out_file)

    # Clean exit
    sys.exit(0)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()