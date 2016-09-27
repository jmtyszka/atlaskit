#!/usr/bin/env python3
"""
Warp a probabistic atlas onto a brain containing a lesion.
Requires a T1w image of the lesioned brain and a binary lesion mask, both in the same space.
This command outputs lesion volumes within each probabilistic atlas label

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
    parser.add_argument('-i', '--input', help='Skull-stripped T1w image of lesioned brain')
    parser.add_argument('-m', '--mask', help='Lesion mask')
    parser.add_argument('-t', '--template', help='Skull-stripped T1w atlas template for registration')
    parser.add_argument('-a', '--atlas', help='Probabilistic atlas labels')

    # Parse command line arguments
    args = parser.parse_args()
    lesion_T1w = args.input
    lesion_mask = args.mask
    atlas_T1w = args.template
    atlas_prob = args.atlas

    # Setup nipype pipeline
    antsnorm = ants.Registration()
    antsnorm.inputs.output_transform_prefix = "trans"
    antsnorm.inputs.fixed_image = lesion_T1w
    antsnorm.inputs.moving_image = atlas_T1w
    antsnorm.inputs.fixed_image_mask = lesion_mask
    antsnorm.inputs.collapse_output_transforms = True
    antsnorm.inputs.initial_moving_transform_com = True
    antsnorm.inputs.num_threads = 4
    antsnorm.inputs.output_inverse_warped_image = True
    antsnorm.inputs.output_warped_image = True
    antsnorm.inputs.sigma_units = ['vox'] * 3
    antsnorm.inputs.transforms = ['Rigid', 'Affine', 'SyN']
    antsnorm.inputs.terminal_output = 'stream'
    antsnorm.inputs.winsorize_lower_quantile = 0.005
    antsnorm.inputs.winsorize_upper_quantile = 0.995
    antsnorm.inputs.convergence_threshold = [1e-06]
    antsnorm.inputs.convergence_window_size = [10]
    antsnorm.inputs.metric = ['MI', 'MI', 'CC']
    antsnorm.inputs.metric_weight = [1.0] * 3
    antsnorm.inputs.number_of_iterations = [[1000, 500, 250, 100], [1000, 500, 250, 100], [100, 70, 50, 20]]
    antsnorm.inputs.radius_or_number_of_bins = [32, 32, 4]
    antsnorm.inputs.sampling_percentage = [0.25, 0.25, 1]
    antsnorm.inputs.sampling_strategy = ['Regular', 'Regular', 'None']
    antsnorm.inputs.shrink_factors = [[8, 4, 2, 1]] * 3
    antsnorm.inputs.smoothing_sigmas = [[3, 2, 1, 0]] * 3
    antsnorm.inputs.transform_parameters = [(0.1,), (0.1,), (0.1, 3.0, 0.0)]
    antsnorm.inputs.use_histogram_matching = True
    antsnorm.inputs.write_composite_transform = True
    antsnorm.run()

    trans = ['trans0GenericAffine.mat', 'trans1InverseWarp.nii.gz']

    apply2anat = ants.ApplyTransforms()
    apply2anat.inputs.input_image = atlas_prob
    apply2anat.inputs.output_image = 'lesion_atlas_prob.nii.gz'
    apply2anat.inputs.transforms = trans
    apply2anat.inputs.default_value = 0
    apply2anat.inputs.dimension = 3
    apply2anat.inputs.input_image_type = 0
    apply2anat.inputs.invert_transform_flags = [True, False]
    apply2anat.inputs.num_threads = 4
    apply2anat.inputs.reference_image = lesion_T1w
    apply2anat.inputs.terminal_output = 'stream'
    apply2anat.run()

    # Clean exit
    sys.exit(0)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
