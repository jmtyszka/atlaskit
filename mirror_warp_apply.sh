#!/bin/bash
# Create bilateral from unilateral labels using precalculated mirror warp
#
# AUTHOR : Mike Tyszka
# PLACE  : Caltech
# DATES  : 2015-07-30 JMT From scratch
#          2016-03-25 JMT Separate calc and apply components of mirror warp
#          2016-03-31 JMT Update for new directory structure

if [ $# -lt 3 ]
then
  echo "USAGE : mirror_warp_apply.sh <Unilateral Labels> <Mirror Warp Directory> <Bilateral Directory>"
  exit
fi

uni_labels=$1
warp_dir=$2
bi_dir=$3

if [ ! -d $warp_dir ]
then
  echo "Mirror warp directory does not exist - exiting"
  exit
fi

# Filenames for unilateral and bilateral label volumes
uni_file=`basename $uni_labels`
uni_stub=${uni_file%%.nii.gz}
bi_labels=${bi_dir}/${uni_file}

# Mirrored version of the unilateral labels (not bilateral yet)
uni_mirror=${bi_dir}/${uni_stub}_mirror.nii.gz
uni_mirror_warp=${bi_dir}/${uni_stub}_mirror_warp.nii.gz

# Bilateral version of labels (original + mirror warped)
labels_bilateral=${bi_dir}/${uni_stub}_bilateral.nii.gz

# ANTs warp prefix
ANTS_prefix=${warp_dir}/Mirror
ANTS_warp=WarpImageMultiTransform
mirror_affine=${ANTS_prefix}Affine.txt
mirror_warp=${ANTS_prefix}Warp.nii.gz

if [ ! -s ${mirror_warp} ]
then
  echo "Mirror warp files do not exist - exiting"
  exit
fi

# Mirrored label image
echo "  Reflecting label image about the mid-sagittal plane"
mirror.py -i ${uni_labels} -o ${uni_mirror}

# Warp flipped LH atlas = unwarped RH atlas to generate correctly warped RH atlas
# It's essential that the original label image is used as a reference here to keep the sform sense correct
# Use nearest neighbour resampling for the label image
echo "  Apply warp to reflected labels"
WarpImageMultiTransform 3 ${uni_mirror} ${uni_mirror_warp} -R ${uni_labels} --use-NN ${mirror_warp} ${mirror_affine}

# Combine LH and RH prob atlases
echo "  Combine original and mirrored warped labels"
fslmaths ${uni_labels} -add ${uni_mirror_warp} ${bi_labels}

# Remove intermediate files
rm -rf ${uni_mirror} ${uni_mirror_warp}
