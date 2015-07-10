#!/bin/bash
#
# Copy left hemisphere labels to the right hemisphere using a diffeomorphic mapping
#
# - Requires ANTs diffeomorphism (affine + warp) mapping LR flipped template to original
# - All labels from the right hemisphere before copying
# - Assumes ANTs used a fixed original atlas template and moving LR flipped template
#
# USAGE : label_copy_lr.sh <atlas> <affine text> <warp image>
#
# OUTPUT : <atlas>_LH.nii.gz, <atlas>_RH.nii.gz and <atlas>_BH.nii.gz
#
# AUTHOR : Mike Tyszka, Ph.D.
# PLACE  : Caltech
# DATES  : 2015-05-05 JMT From scratch

if [ $# -lt 4 ]; then
  echo "USAGE : label_copy_lr.sh <source atlas> <destination atlas image> <affine text> <warp image>"
  exit
fi

# Parse arguments
atlas=$1
affine=$2
warp=$3

# Output filenames
stub=${atlas%%.nii.gz}
atlas_LH=${stub}_LH.nii.gz
atlas_RH=${stub}_RH.nii.gz
atlas_BH=${stub}_BH.nii.gz

# Clear right hemisphere (assumes R-L in x)
hx=`fslinfo $src_atlas | awk '{ if ($1 == "dim1") print $2/2 }'`
echo "  Clearing atlas right hemisphere"
fslroi $atlas $atlas_LH 0 $hx 0 -1 0 -1 0 -1

# Flip new atlas L-R
echo "  Flipping atlas L-R"
fslswapdim $atlas_LH -x y z $atlas_RH

# Apply diffeomorphism to new atlas. NN interpolation
echo "  Warping atlas right hemisphere"
WarpImageMultiTransform 3 ${atlas_RH} ${atlas_RH} -R ${atlas} --use-NN $warp $affine

# Add LH and RH atlases to create BH atlas
echo "  Combining left and right hemispheres"
fslmaths $atlas_LH -add $atlas_RH $atlas_BH

echo "Results in $atlas_LH, $atlas_RH and $atlas_BH"




