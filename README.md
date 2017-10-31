# atlaskit
Tools for working with anatomic atlas labels.

## Introduction
These python and bash scripts provide convenience features for working with integer-valued and probabilistic label images for digital human brain atlases. Some of the scripts wrap a combination of FSL and ANTs commands which map individual T1w and/or T2w structural images to most high SNR midspace templates, including MNI/ICBM152 and CIT168/HCP.

## Requirements
| Software | Description | Version | Link |
| :------- | :---------- | :------ | :--- |
| Python   | High-level coding | 3.4+ | https://www.python.org/ |
| FSL      | Neuroimaging analysis | 5.0.9+ | https://fsl.fmrib.ox.ac.uk/ |
| ANTs     | Image warp registration | 2.1.0+ | https://github.com/stnava/ANTs/ |
| CIT168   | Templates and amygdala atlas | 1.0.1+ | http://evendim.caltech.edu/amygdala-atlas/ |

## Registering the CIT168 atlas to an individual brain
The CIT168 atlas provides both a T1w and T2w template. Both are accurately registerd to each other (ie they're in the same space) so you can register the atlas to just an individual T1w image or to a similar pair of T1w and T2w images from an individual. If you only have T1w 3D structural images, use the tmp2ind_T1.sh script. If you have both T1w and T2w 3D structural images for an individual participant, use the tmp2ind_T1T2.sh script which makes use of both contrasts during the registration.

Run either script with no arguments for help on usage.

## Creating a midspace template from your own data
We've included a pair of scripts for convenience that simply wrap the antsMultivariateTemplateConstruction.sh command. midspace_T1.sh constructs a minimum deformation midspace template from a set of T1w individual images only, and midspace_T1T2.sh constructs a T1w and T2w template pair from accurately registered T1w and T2w individual images.

To run either script:

1. Create a directory containing all your individual images.

2. If you have both T1w and T2w images for each individual, make sure that

  (a) They are accurately registered and resampled to the same resolution and field of view. This can be achieved using any linear registration command that supports a mutual information metric (eg FLIRT).

  (b) The subject ID appears before the contrast type (T1 or T2) in the filename. For example the image pair for subject 123456 could be named sub-123456_T1w_brain.nii.gz and sub-123456_T2w_brain.nii.gz.

3. Copy the required script to the image directory. Make sure it is executable using chmod +x *.sh and run from the image directory.

Template generation could take from several hours to several days depending on the number of individual images, image dimensions and computing resources available.

## Bugs, Issues, Feature Requests

If you encounter any problems using these scripts or the CIT168 atlas, let us know by raising a new issue for this repository through https://github.com/jmtyszka/atlaskit/issues
