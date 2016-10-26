#!/usr/bin/env python3
"""
Warp a probabistic atlas onto a brain containing a lesion.
Requires a T1w image of the lesioned brain and a corresponding binary lesion mask in the same space.

Usage
----
atlas2lesion.py -i <lesion T1w> -m <lesion mask> -t <atlas T1w template> -a <probabilistic atlas>
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
from nipype.interfaces import ants


def main():

    # Construct a command line argument parser
    parser = argparse.ArgumentParser(description='Atlas-based lesion volumetrics')
    parser.add_argument('-i', '--input', required=True, help='Skull-stripped T1w image of lesioned brain')
    parser.add_argument('-m', '--mask', required=True, help='Lesion mask')
    parser.add_argument('-t', '--template', required=True, help='Skull-stripped T1w atlas template for registration')
    parser.add_argument('-a', '--atlas', required=True, help='Probabilistic atlas labels')

    # Parse command line arguments
    args = parser.parse_args()
    lesion_T1w = args.input
    lesion_mask = args.mask
    atlas_T1w = args.template
    atlas_prob = args.atlas

    # ANTS composite (affine + SyN warp) transform file
    xfms = ['trans0GenericAffine.mat','trans1Warp.nii.gz']

    # Run ANTs SyN registration of atlas T1w template to T1w lesion image
    if not os.path.isfile(xfms[0]):
        print('  Calculating warp from template to individual space')
        antsreg = ants.Registration()
        antsreg.inputs.output_transform_prefix = "trans"
        antsreg.inputs.fixed_image = lesion_T1w
        antsreg.inputs.moving_image = atlas_T1w
        antsreg.inputs.fixed_image_mask = lesion_mask
        antsreg.inputs.collapse_output_transforms = True
        antsreg.inputs.initial_moving_transform_com = True
        antsreg.inputs.num_threads = 4
        antsreg.inputs.smoothing_sigmas = [[3, 2, 1, 0]] * 3
        antsreg.inputs.sigma_units = ['vox'] * 3
        antsreg.inputs.transforms = ['Rigid', 'Affine', 'SyN']
        antsreg.inputs.terminal_output = 'stream'
        antsreg.inputs.winsorize_lower_quantile = 0.005
        antsreg.inputs.winsorize_upper_quantile = 0.995
        antsreg.inputs.convergence_threshold = [1e-06]
        antsreg.inputs.convergence_window_size = [10]
        antsreg.inputs.metric = ['MI', 'MI', 'CC']
        antsreg.inputs.metric_weight = [1.0] * 3
        antsreg.inputs.number_of_iterations = [[1000, 500, 250, 100], [1000, 500, 250, 100], [100, 70, 50, 20]]
        antsreg.inputs.radius_or_number_of_bins = [32, 32, 4]
        antsreg.inputs.sampling_strategy = ['Regular', 'Regular', 'None']
        antsreg.inputs.sampling_percentage = [0.25, 0.25, 1]
        antsreg.inputs.shrink_factors = [[8, 4, 2, 1]] * 3
        antsreg.inputs.transform_parameters = [(0.1,), (0.1,), (0.1, 3.0, 0.0)]
        antsreg.inputs.use_histogram_matching = True
        antsreg.inputs.write_composite_transform = False
        antsreg.inputs.output_warped_image = True
        antsreg.inputs.output_inverse_warped_image = True
        antsreg.run()
    else:
        print('* ANTs registration already run - continuing to resampling')

    # Note that the affine+SyN warp is written as a single composite transform
    if os.path.isfile(xfms[0]):
        print('  Warping probabilistic atlas to individual space')
        warp = ants.WarpTimeSeriesImageMultiTransform()
        warp.inputs.input_image = atlas_prob
        warp.inputs.out_postfix = '_warped'
        warp.inputs.transformation_series = xfms
        warp.inputs.num_threads = 4
        warp.inputs.reference_image = lesion_T1w
        warp.inputs.terminal_output = 'stream'
        warp.run()
    else:
        print('* ANTs composite transform file not found - exiting')

    # Clean exit
    sys.exit(0)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
