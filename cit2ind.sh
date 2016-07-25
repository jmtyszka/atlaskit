#!/bin/bash
# Register CIT168 atlas to an individual T1w image
#
# REQUIRES:
# 1. Working installation of the ANTs registration package, specifically:
# - N4BiasFieldCorrection
# - antsApplyTransforms
# - antsRegistrationSyNQuick.sh (from the Scripts release)
# 2. Working installation of FSL 5.0.8
#
# AUTHOR : Mike Tyszka
# PLACE  : Caltech
# DATES  : 2015-12-09 JMT Adapt ANTs brain extraction script to use exact brain template
#          2015-12-16 JMT Adapt brainiac to extract brain, then warp CIT168 atlas to individual space
#
# Copyright 2015 California Institute of Technology

if [ $# -lt 1 ]
then
	echo "USAGE : cit2ind <T1 whole-head image>"
	exit
fi

# Get T1 image filename
T1=$1

if [ ! -e $T1 ]
then
	echo "*** ${T1} does not exist - exiting"
	exit
fi

# Check for CIT_ATLAS_DIR env variable
if [ -z "$CIT_ATLAS_DIR" ]
then
	echo "Environment variable CIT_ATLAS_DIR needs to be defined"
	exit
fi

# MNI 3 mm head template from CIT Atlas release
T1_MNI_3mm=${CIT_ATLAS_DIR}/MNI152/MNI152_T1_3mm.nii.gz
T1_MNI_3mm_mask=${CIT_ATLAS_DIR}/MNI152/MNI152_T1_3mm_CIT_brain_mask.nii.gz

echo ""
echo "BRAINIAC"
echo "--------"
echo "Extracting brain from ${T1}"

# Create stub without .nii or .nii.gz extension
stub=${T1%%.gz}
stub=${stub%%.nii}
echo "  Output stub is ${stub}"

# Scratch directory
scratch_dir=${stub}.brainiac
mkdir -p ${scratch_dir}

# Intermediate filenames
T1_3mm=${scratch_dir}/${stub}_3mm.nii.gz
T1_3mm_nobias=${scratch_dir}/${stub}_3mm_nobias.nii.gz
T1_3mm_bias=${scratch_dir}/${stub}_3mm_bias.nii.gz
T1_brain=${scratch_dir}/${stub}_brain.nii.gz
T1_brain_mask=${scratch_dir}/${stub}_brain_mask.nii.gz

# Warp filenames
Prefix=${scratch_dir}/MNI2Ind
Affine=${Prefix}0GenericAffine.mat
Warp=${Prefix}1Warp.nii.gz

# Downsample T1 image to 3 mm isotropic
echo "  Downsampling to 3 mm"
if [ ! -e ${T1_3mm} ]
then
	flirt -in ${T1} -out ${T1_3mm} -ref ${T1} -init ${FSLDIR}/etc/flirtsch/ident.mat -applyisoxfm 3.0 -interp sinc
fi

# N4 bias correct T1 image
echo "  N4 bias correcting"
if [ ! -e ${T1_3mm_nobias} ]
then
	N4BiasFieldCorrection -d 3 -i ${T1_3mm} -o [${T1_3mm_nobias},${T1_3mm_bias}]
fi

# Warp to 3 mm MNI head template
echo "  Warping MNI to individual space"
if [ ! -e ${Warp} ]
then
	antsRegistrationSyNQuick.sh -d 3 -f ${T1_3mm_nobias} -m ${T1_MNI_3mm} -o ${Prefix} > /dev/null
fi

# Apply warp to 3 mm brain mask and resample to original T1 space
echo "  Warping brain mask to individual space"
if [ ! -e ${T1_brain_mask} ]
then
	antsApplyTransforms -d 3 -i ${T1_MNI_3mm_mask} -r ${T1} -o ${T1_brain_mask} -n BSpline -t ${Warp} -t ${Affine}
fi

# Create brain-only image by multiplying the original T1 and the soft mask
echo "  Soft masking brain"
if [ ! -e ${T1_brain} ]
then
	fslmaths ${T1} -mul ${T1_brain_mask} ${T1_brain}
fi

echo "Complete - all results are in ${scratch_dir}"
