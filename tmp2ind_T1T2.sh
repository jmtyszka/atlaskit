#!/bin/bash
#
# Perform SyN multimodal warping from individual to template space using joint T1 and T2 cost function
#
# AUTHOR : Mike Tyszka
# PLACE  : Caltech
# DATES  : 2016-09-30 JMT From scratch
#          2017-04-10 JMT Fixed dimensions bug in pAtlas resampling
#          2017-04-11 JMT Duplicated Adam Mezher's fixes from T1 to T1T2 version
#          2017-06-13 JMT Simplied argument handling
#	   2018-04-20 JD  added names for arguments; added nthreads as argument;
#                         using antsApplyTransforms instead of WarpImageMultiTransform; made atlas argument optional
#
# MIT License
#
# Copyright (c) 2017 Caltech
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

# function for parsing options
getopt1() {

    sopt="$1"
    shift 1
    for fn in $@ ; do
    if [ `echo $fn | grep -- "^${sopt}=" | wc -w` -gt 0 ] ; then
        echo $fn | sed "s/^${sopt}=//"
        return 0
    fi
    done
}

if [ $# -lt 4 ]
then
    echo "USAGE : tmp2ind_T1T2.sh" 
    echo "        [REQUIRED ARGS] --T1ind=<Individual T1w> --T2ind=<Individual T2w> --T1tmp=<Template T1w image> --T2tmp=<Template T2w image>"
    echo "        [OPTIONAL ARGS] --Atmp=<(Probabilistic) atlas=''> --outdir=<output directory> --nthreads=<nthreads=4>"
    echo "All images are Nifti-1 format, compressed or uncompressed"
    echo "This script should be called from within the directory that you wish transforms to be saved to"
    exit
fi


T1ind=`getopt1 "--T1ind" $@`
if [ ! -s $T1ind ]; then
    echo "*T1ind: $T1ind does not exist or is empty"
    exit
fi

T2ind=`getopt1 "--T2ind" $@`
if [ ! -s $T2ind ]; then
    echo "*T2ind: $T2ind does not exist or is empty"
    exit
fi

T1tmp=`getopt1 "--T1tmp" $@`
if [ ! -s $T1tmp ]; then
    echo "*T1tmp: $T1tmp does not exist or is empty"
    exit
fi

T2tmp=`getopt1 "--T2tmp" $@`
if [ ! -s $T2tmp ]; then
    echo "*T2tmp: $T2tmp does not exist or is empty"
    exit
fi

pAtmp=`getopt1 "--Atmp" $@`
if [ -z "$pAtmp" ]; then
    echo "*Atmp is empty"
elif [ ! -s $pAtmp ]; then
    echo "*Atmp: $pAtmp does not exist"
fi

outdir=`getopt1 "--outdir" $@`
if [ -z "$outdir" ];then
    echo "*outdir is empty"
    outdir=TMP2IND
fi

nthreads=`getopt1 "--nthreads" $@`
if [ -z "$nthreads" ];then
    echo "*nthreads is empty; set to default value = 4"
    nthreads=4
fi

# Splash text
echo ""
echo "------------------------------------------------------------"
echo " SyN Warp T1 and T2 templates to individual space"
echo "------------------------------------------------------------"
echo "Individual T1 : ${T1ind}"
echo "Individual T2 : ${T2ind}"
echo "  Template T1 : ${T1tmp}"
echo "  Template T2 : ${T2tmp}"
if [ ! -z "$pAtmp" ] && [ -s ${pAtmp} ]
then
    echo "   Prob Atlas : ${pAtmp}"
fi
echo "Output directory : ${outdir}"
echo "Using ${nthreads} threads"

# Create output directory
mkdir -p ${outdir}

# Registration files
prefix=${outdir}/TMP2IND
tmp2ind_affine=${prefix}0GenericAffine.mat
tmp2ind_warp=${prefix}1Warp.nii.gz
logfile=${prefix}_Warp.log

# Output filenames
T1tmp2ind=${outdir}/T1w_tmp2ind.nii.gz
T2tmp2ind=${outdir}/T2w_tmp2ind.nii.gz
pAtmp2ind=${outdir}/pA_tmp2ind.nii.gz
dAtmp2ind=${outdir}/dA_tmp2ind.nii.gz

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
	antsApplyTransforms -d 3 -e 0 -i ${T2tmp} -r ${T1ind} -o ${T2tmp2ind} -n BSpline -t ${tmp2ind_warp} -t ${tmp2ind_affine}
fi

if [ ! -z "$pAtmp" ] && [ -s ${pAtmp} ] && [ ! -s ${pAtmp2ind} ]
then
    echo "Warping probabilistic atlas into individual space"
    antsApplyTransforms -d 3 -e 3 -i ${pAtmp} -r ${T1ind} -o ${pAtmp2ind} -n NearestNeighbor -t ${tmp2ind_warp} -t ${tmp2ind_affine}
fi

if [ ! -s ${dAtmp2ind} ]
then
    echo "Creating deterministic labels in individual space (p > 0.25)"
    prob2det.sh ${pAtmp2ind} ${dAtmp2ind} 0.25
fi

# Report output filenames
echo ""
echo "------------------------------------------------------------"
echo " Output files"
echo "------------------------------------------------------------"
echo "Template T1 in individual space : ${T1tmp2ind}"
echo "Template T2 in individual space : ${T2tmp2ind}"

if [ ! -z "$pAtmp" ] && [ -s ${pAtmp} ] && [ -s ${pAtmp2ind} ]
then
    echo "Prob atlas in individual space  : ${pAtmp2ind}"
fi

if [ -s ${dAtmp2ind} ]
then
    echo "Det atlas in individual space   : ${dAtmp2ind}"
fi
