#!/bin/bash
# Register CIT168 atlas to MNI image
#
# REQUIRES:
# 1. Working installation of the ANTs registration package
#
# AUTHOR : Wolfgang M. Pauli 
# PLACE  : Caltech
# DATES  : 2017-08-20 WMP, from scratch
#
# Copyright 2017 California Institute of Technology


if [ $# -lt 1 ]
then
	echo "USAGE : cit2mni.sh <MNI resolution in mm (1 or [2])> <n threads>"
	exit
fi

# set desired resolution
res=$1

# set desired n of threads
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$2

# Check for FSL installation
if [ -z  $FSLDIR ]
then
	echo "Environment variable FSLDIR needs to be defined"
	exit
fi

# Check for ANTS installation
if [ -z  $ANTSPATH ]
then
	echo "Environment variable ANTSPATH needs to be defined"
	exit
fi


# Check for CIT_ATLAS_DIR env variable
if [ -z "$CIT_ATLAS_DIR" ]
then
	echo "Environment variable CIT_ATLAS_DIR needs to be defined"
	exit
fi


T1_MNI=${FSLDIR}/data/standard/MNI152_T1_${res}mm_brain.nii.gz
T1_CIT=${CIT_ATLAS_DIR}/CIT168_700um/CIT168_T1w_700um.nii.gz

# Create brain-only image by multiplying the original T1 and the soft mask
if [ ! -e ${T1_MNI} ]
then
    echo "Could not find fixed image ${T1_MNI}"
else
    echo "Fixed image: ${T1_MNI}"
fi

if [ ! -e ${T1_CIT} ]
then
    echo "Could not find moving image ${T1_CIT}"
else
    echo "Moving image: ${T1_CIT}"
fi


# Mattes metric parameters
metricWeight=1
radius=4
dim=3

# echo  "Running affine initializer"
# ${ANTSPATH}antsAffineInitializer $dim $T1_CIT $T1_MNI ${CIT_ATLAS_DIR}/CIT_2_MNI_${res}mm_init.mat

# exit

echo  "Running ants registration"
# Deformable
${ANTSPATH}/antsRegistration \
   --verbose 1 \
   --dimensionality $dim \
   --metric CC[ $T1_MNI, $T1_CIT, $metricWeight, $radius] \
   --transform Syn[0.1,3,0] \
   --convergence [100x100x100x100x100, 1.e-6, 10] \
   --smoothing-sigmas 5x3x2x1x0vox \
   --shrink-factors 10x6x4x2x1 \
   --interpolation Linear \
   -o [ ${CIT_ATLAS_DIR}/CIT168toMNI152_${res}mm/CIT168_${res}mm_MNI, ${CIT_ATLAS_DIR}/CIT168toMNI152_${res}mm/CIT168_${res}mm_MNI_warped.nii.gz ]

#   --use-histogram-matching 1 \
# -x [ $fixed_mask, $moving_mask ] 
#    --initial-moving-transform ${CIT_ATLAS_DIR}/CIT_2_MNI_${res}mm_init.mat \   
