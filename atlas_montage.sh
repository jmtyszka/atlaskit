#!/bin/bash
# Create overlay montage of a probabilistic atlas on tissue boundaries from the associate structural template
#
# AUTHOR : Mike Tyszka
# PLACE  : Caltech
# DATES  : 2017-05-23 JMT From scratch

if [ $# -lt 3 ]
then
  echo "USAGE : atlas_montage.sh <4D prob atlas> <3D structural template> <p threshold>"
  exit
fi

PROB=$1
STRUCT=$2
THRESH=$3

# Colormap for deterministic atlas overlay
LUT=${FSLDIR}/etc/luts/striatum-con-7sub.lut

# Temporary filenames
TMP_MASK=`mktemp`_mask.nii.gz
TMP_MASK_CROP=`mktemp`_mask_crop.nii.gz
TMP_PROB_CROP=`mktemp`_prob_crop.nii.gz
TMP_DET_CROP=`mktemp`_det_crop.nii.gz
TMP_STRUCT_CROP=`mktemp`_struct_crop.nii.gz

# Create binary mask for suprathreshold voxels over all labels
echo "Creating mask for suprathreshold voxels over all labels"
fslmaths ${PROB} -Tmax -thr ${THRESH} ${TMP_MASK}

# Find minimum bounding box
BB_3D=`fslstats ${TMP_MASK} -w`
BB_4D=${BB_3D/0 1/0 -1}
echo "3D bounding box set to ${BB_3D}"
echo "4D bounding box set to ${BB_4D}"

# Apply bounding box to atlas, structural and mask
echo "Cropping probabilistic atlas"
fslroi ${PROB} ${TMP_PROB_CROP} ${BB_4D}
echo "Cropping structural template"
fslroi ${STRUCT} ${TMP_STRUCT_CROP} ${BB_3D}
echo "Cropping mask"
fslroi ${TMP_MASK} ${TMP_MASK_CROP} ${BB_3D}

# Create cropped deterministic atlas
echo "Creating deterministic atlas"
fslmaths ${TMP_PROB_CROP} -Tmaxn -add 1 -mas ${TMP_MASK_CROP} ${TMP_DET_CROP}

# Flip dimensions to make slicer output coronal instead of axial
echo "Reorienting dimensions for coronal sections"
fslswapdim ${TMP_DET_CROP} x z -y ${TMP_DET_CROP}
fslswapdim ${TMP_STRUCT_CROP} x z -y ${TMP_STRUCT_CROP}

# Overlay deterministic atlas on edges from structural
echo "Overlaying deterministic atlas on structural reference"
slicer ${TMP_DET_CROP} ${TMP_STRUCT_CROP} -l ${LUT} -S 4 256 montage.png

# Cleanup
echo "Cleaning up"
rm -rf ${TMP_MASK} ${TMP_MASK_CROP} ${TMP_PROB_CROP} ${TMP_DET_CROP} ${TMP_STRUCT_CROP}
