#!/bin/bash
# Convert probabilistic 4D atlas (labels in 4th dimension) to indexed 3D volume with thresholding
# - Use approach favored by FSL
#
# AUTHOR : Mike Tyszka
# PLACE  : Caltech
# DATES  : 2015-08-19 JMT From scratch
#          2016-04-20 JMT Use mktemp filename for mask to avoid SGE overwrites

if [ $# -lt 3 ]
then
  echo "USAGE : prob2ind.sh <4D prob atlas> <3D label atlas> <p threshold>"
  exit
fi

PROB_ATLAS=$1
IND_ATLAS=$2
THRESH=$3

# Temporary mask for this prob atlas
MASK=${PROB_ATLAS%%.nii.gz}_tmp.nii.gz

fslmaths ${PROB_ATLAS} -Tmax -thr ${THRESH} ${MASK}
fslmaths ${PROB_ATLAS} -Tmaxn -add 1 -mas ${MASK} ${IND_ATLAS}

# Cleanup
rm -rf ${MASK}
