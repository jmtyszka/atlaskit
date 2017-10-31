#!/bin/bash
#
# Perform SyN multimodal warping from individual to template space using joint T1 and T2 cost function
#
# USAGE  : tmp2ind_T1T2.sh <Individual T1w> <Individual T2w> <Template T1w image> <Template T2w image> <Prob atlas>"
#
# AUTHOR : Mike Tyszka
# PLACE  : Caltech
# DATES  : 2016-09-30 JMT From scratch
#          2017-04-10 JMT Fixed dimensions bug in pAtlas resampling
#          2017-04-11 JMT Duplicated Adam Meher's fixes from T1 to T1T2 version
#          2017-06-13 JMT Simplied argument handling
#
# MIT License
#
# Copyright (c) 2017 Mike Tyszka
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

if [ $# -lt 5 ]
then
	echo "USAGE : tmp2ind_T1T2.sh <Individual T1w> <Individual T2w> <Template T1w image> <Template T2w image> <Probabilistic atlas>"
	echo "All images are Nifti-1 format, compressed or uncompressed"
	exit
fi


T1ind=$1
if [ ! -s $T1ind ]; then
    echo "* $T1ind does not exist or is empty"
    exit
fi

T2ind=$2
if [ ! -s $T2ind ]; then
    echo "* $T2ind does not exist or is empty"
    exit
fi

T1tmp=$3
if [ ! -s $T1tmp ]; then
    echo "* $T1tmp does not exist or is empty"
    exit
fi

T2tmp=$4
if [ ! -s $T2tmp ]; then
    echo "* $T2tmp does not exist or is empty"
    exit
fi

pAtmp=$5
if [ ! -s $pAtmp ]; then
    echo "* $pAtmp does not exist or is empty"
    exit
fi

# Splash text
echo "------------------------------------------------------------"
echo " SyN Warp T1 and T2 templates to individual space"
echo "------------------------------------------------------------"
echo "Individual T1 : ${T1ind}"
echo "Individual T2 : ${T2ind}"
echo "  Template T1 : ${T1tmp}"
echo "  Template T2 : ${T2tmp}"
echo "   Prob Atlas : ${pAtmp}"

# Fixed ANTs parameters
nthreads=4

# Registration files
prefix=TMP2IND_
tmp2ind_affine=${prefix}0GenericAffine.mat
tmp2ind_warp=${prefix}1Warp.nii.gz
logfile=${prefix}Warp.log

# Output filenames
T1tmp2ind=T1w_tmp2ind.nii.gz
T2tmp2ind=T2w_tmp2ind.nii.gz
pAtmp2ind=pA_tmp2ind.nii.gz

# Calculate affine and SyN warp
if [ ! -s ${tmp2ind_warp} ]
then
    echo "Starting SyN registration"
	antsRegistrationSyN.sh -d 3 -n ${nthreads} -t b -o ${prefix} -f ${T1ind} -f ${T2ind} -m ${T1tmp} -m ${T2tmp} 2>&1 > ${logfile}
fi

# Rename warped template T1
if [ ! -s ${T1tmp2ind} ]
then
	mv ${prefix}Warped.nii.gz ${T1tmp2ind}
fi

# Resample template T2 to individual space
if [ ! -s ${T2tmp2ind} ]
then
    echo "Warping T2w template into individual space"
	WarpImageMultiTransform	3 ${T2tmp} ${T2tmp2ind} -R ${T1ind} ${tmp2ind_warp} ${tmp2ind_affine} --use-BSpline
fi

# Resample probabilistic atlas to individual space
if [ ! -s ${pAtmp2ind} ]
then
    echo "Warping probabilistic atlas into individual space"
	WarpImageMultiTransform	4 ${pAtmp} ${pAtmp2ind} -R ${T1ind} ${tmp2ind_warp} ${tmp2ind_affine} --use-BSpline
fi

# Report output filenames
echo "------------------------------------------------------------"
echo " Output files"
echo "------------------------------------------------------------"
echo "Template T1 in individual space : ${T1tmp2ind}"
echo "Template T2 in individual space : ${T2tmp2ind}"
echo "Prob atlas in individual space  : ${pAtmp2ind}"
