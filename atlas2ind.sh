#!/bin/bash
#
# Warp atlas labels from template (350um) to individual spaces
# Requires previous ANTs diffeomorphic registration of individual to atlas space (or vise versa)
#
# USAGE: atlas2ind.sh <Affine> <Warp> <Labels> <Output> <Direction>
#
# ARGS:
# Affine    : Affine text file (ANTs output) 
# Warp      : Warp field image. Use Warp for dirn = 0 and InverseWarp for dirn = 1
# Labels    : Atlas labels image
# Output    : Output label file name (in individual space)
# Direction : Original mapping direction generating affine and warp fields
#             0 = atlas to individual, 1 = individual to atlas
#
# AUTHOR : Mike Tyszka
# PLACE  : Caltech
# DATES  : 2015-04-28 JMT Adapt from indlabel.sh

if [ $# -lt 5 ]; then
  echo "USAGE: atlas2ind.sh <Affine> <Warp> <Labels> <Output> <Direction>"
  exit
fi

# Assign arguments
Affine=$1
Warp=$2
Labels=$3
Output=$4
Direction=$5

# Warp atlas to individual space with nearest neighbor interpolation
if [ $Direction -eq "0" ]; then
  # Original was atlas to individual mapping. No need to invert.
  WarpImageMultiTransform 3 ${ATLAS} ${ind_atlas} -R ${img} --use-NN ${Warp} ${Affine} > /dev/null 2>&1
elif [ $Direction -eq "1" ]; then
  # Original was individual to atlas mapping. Invert order of affine and warp 
  WarpImageMultiTransform 3 ${ATLAS} ${ind_atlas} -R ${img} --use-NN -i ${Affine} ${Warp} > /dev/null 2>&1
else
  echo "Unknown mapping direction: type atlas2ind.sh for usage."
  exit
fi
