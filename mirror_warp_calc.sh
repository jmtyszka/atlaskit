#!/bin/bash
# Calculate the diffeomorphic warp required to map one hemisphere onto the other.
# Requires an initial reflection of images about the midsagittal plane followed by the
# diffeomorphic warping of the reflected to the unreflected image.
#
# This implementation requires both T1 and T2 contrasts for a multivariate optimization.
#
# AUTHOR : Mike Tyszka
# PLACE  : Caltech
# DATES  : 2015-07-30 JMT From scratch
#          2016-03-25 JMT Separate calc and apply components of mirror warp

if [ $# -lt 3 ]
then
  echo "USAGE : mirror_warp_calc.sh <T1> <T2> <Output Directory>"
  exit
fi

T1_orig=$1
T2_orig=$2
out_dir=$3

# Create output directory safely
mkdir -p $out_dir

# Filename stubs for original images
T1_fname=`basename $T1_orig`
T1_stub=${T1_fname%%.nii.gz}
T2_fname=`basename $T2_orig`
T2_stub=${T2_fname%%.nii.gz}

# Copies of original T1 and T2 in output directory
T1=${out_dir}/${T1_fname}
T2=${out_dir}/${T2_fname}

# LR mirrored versions of original images
T1_mirror=${out_dir}/${T1_stub}_mirror.nii.gz
T2_mirror=${out_dir}/${T2_stub}_mirror.nii.gz

# Final mirror warped versions of original images
T1_mirror_warp=${out_dir}/${T1_stub}_mirror_warp.nii.gz
T2_mirror_warp=${out_dir}/${T2_stub}_mirror_warp.nii.gz

# ANTs output prefix
ANTS_prefix=${out_dir}/Mirror
ANTS_warp=WarpImageMultiTransform
mirror_affine=${ANTS_prefix}Affine.txt
mirror_warp=${ANTS_prefix}Warp.nii.gz

# Bivariate ANTs registration options
ANTS_OPTS="-i 30x90x20 -t SyN[0.25] -r Gauss[3,0] --use-Histogram-Matching --number-of-affine-iterations 10000x10000x1000 --MI-option 32x16000"

# Copy T1 and T2 images into output directory
echo "  Copying original images to output directory"
cp $T1_orig $T1
cp $T2_orig $T2
chmod +w $T1 $T2

# L-R flip master T1 and T2 templates
echo "  Reflecting images about the mid-sagittal plane"
mirror.py -i $T1 -o $T1_mirror
mirror.py -i $T2 -o $T2_mirror

if [ ! -s ${mirror_warp} ]
then

  # Warp L-R flipped T1 and T2 onto their unflipped versions
  echo "  Warping L-R reflected images onto their unreflected versions"
  ANTS 3 -m CC[${T1},${T1_mirror},1,5] -m CC[${T2},${T2_mirror},1,5] ${ANTS_OPTS} -o ${ANTS_prefix}

	# Apply mirror warp to L-R flipped T1 and T2
	echo "  Applying mirror warp"
	WarpImageMultiTransform 3 ${T1_mirror} ${T1_mirror_warp} $mirror_warp $mirror_affine
	WarpImageMultiTransform 3 ${T2_mirror} ${T2_mirror_warp} $mirror_warp $mirror_affine

else

  echo "  Previous version of warp exists - exiting"

fi