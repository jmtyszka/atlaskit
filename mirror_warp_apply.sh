#!/bin/bash
# Create bilateral probabilistic labels from unilateral labels using precalculated mirror warp
#
# AUTHOR : Mike Tyszka
# PLACE  : Caltech
# DATES  : 2015-07-30 JMT From scratch
#          2016-03-25 JMT Separate calc and apply components of mirror warp
#          2016-03-31 JMT Update for new directory structure
#          2017-10-29 JMT Use T1 in Warp directory as reference

if [ $# -lt 3 ]
then
  echo "USAGE : mirror_warp_apply.sh <4D Unilateral Probabilistic Labels> <Mirror Warp Directory> <Bilateral Directory>"
  exit
fi

uni_labels=$1
warp_dir=$2
bi_dir=$3

echo "-----------------------------------------"
echo "Apply mirror warp to probabilistic labels"
echo "-----------------------------------------"

if [ ! -d $warp_dir ]
then
  echo "*** Mirror warp directory not found - exiting"
  exit
fi

# Unilateral label filename and stub
uni_file=`basename $uni_labels`
uni_stub=${uni_file%%.nii.gz}

# Mirrored version of the unilateral prob labels (not bilateral yet)
uni_mirror=${bi_dir}/${uni_stub}_mirror.nii.gz
uni_mirror_warp=${bi_dir}/${uni_stub}_mirror_warp.nii.gz

# Bilateral version of prob labels (original + mirror warped)
bi_labels=${bi_dir}/${uni_stub}_bilateral.nii.gz

# ANTs warp prefix
ANTS_prefix=${warp_dir}/Mirror
ANTS_warp=WarpImageMultiTransform
mirror_affine=${ANTS_prefix}Affine.txt
mirror_warp=${ANTS_prefix}Warp.nii.gz
reference=${warp_dir}/T1.nii.gz

if [ ! -s ${mirror_warp} ]
then
  echo "Mirror warp files do not exist - exiting"
  exit
fi

# LR mirrored prob labels
echo "  Reflecting label image about the mid-sagittal plane"
mirror.py -i ${uni_labels} -o ${uni_mirror}

# Warp flipped LH atlas = unwarped RH atlas to generate correctly warped RH atlas
# It's essential that the original 3D T1 template is used as a reference here to keep the sform sense correct
echo "  Apply warp to reflected prob labels"
WarpTimeSeriesImageMultiTransform 4 ${uni_mirror} ${uni_mirror_warp} -R ${reference} ${mirror_warp} ${mirror_affine}

# Combine LH and RH prob atlases
echo "  Combine original and mirrored warped prob labels"
fslmaths ${uni_labels} -add ${uni_mirror_warp} ${bi_labels}

# Remove intermediate files
# rm -rf ${uni_mirror} ${uni_mirror_warp}
