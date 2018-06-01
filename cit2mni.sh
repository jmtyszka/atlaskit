#!/bin/bash
# ANTs B-spline SyN warp CIT168 head template(s) to an MNI space
# - currently supports FSL MNI152 1 mm or MNI ICBM152 asym nonlin 2009c 1 mm spaces
# - uses whole head templates for warping
# - T2w templates used if available (eg with 2009c)
#
# USAGES
# AUTHOR : Mike Tyszka
# PLACE  : Caltech
# DATES  : 2017-10-31 JMT Adapt from mirror_warp_calc.sh
#          2018-05-30 JMT Generalize to other MNI spaces
#
# License
# ----
# This file is part of atlaskit.
#
#     atlaskit is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     atlaskit is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with atlaskit.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright
# ----
# 2018 California Institute of Technology.

# Key directories (edit as needed)
cit_dir="${HOME}/Box/CIT_SubCortical_Atlas/Labeling/reinf_learn_atlas"
icmb_dir="${HOME}/Box/Atlases/MNI_ICBM152_2009c/mni_icbm152_nlin_asym_09c"

# CIT168 700 um templates and prob atlas
cit_t1=${cit_dir}/CIT168_T1w_head_700um.nii.gz
cit_t2=${cit_dir}/CIT168_T2w_head_700um.nii.gz

# Reinforcement learning prob atlas
RL_prob="${HOME}/Box

if [ $# -lt 1 ]
then
    echo "USAGE : cit2mni.sh <atlas space>"
    echo "EXAMPLE : cit2mni.sh 2009c"
    echo ""
    echo "Supported atlas spaces:"
    echo "  FSL   : MNI152 6th generation T1w template used by FSL, SPM"
    echo "  2009c : ICMB152 2009c nonlinear asymmetric T1w and T2w templates"
fi

atlas_space=$1

case ${atlas_space} in
    "FSL")
        out_space="MNI152-FSL"
        mni_t1=${FSLDIR}/data/standard/MNI152_T1_1mm.nii.gz
        mni_t2=""
        ;;
    "2009c")
        out_space="MNI152-2009c"
        mni_t1=${icmb_dir}/mni_icbm152_t1_tal_nlin_asym_09c.nii
        mni_t2=${icmb_dir}/mni_icbm152_t2_tal_nlin_asym_09c.nii
        ;;
    *)
        echo "* Unknown atlas space - exiting"
        exit
esac

# Output directory (subdir of current dir)
out_dir=./CIT168to${out_space}
mkdir -p ${out_dir}

# ANTs output prefix and transform files
ants_prefix=${out_dir}/cit2mni_
cit2mni_affine=${ants_prefix}0GenericAffine.mat
cit2mni_warp=${ants_prefix}1Warp.nii.gz

# ANTs logfile
logfile=${ants_prefix}Warp.log
echo "ANTs registration output in ${logfile}"

# CIT168 templates in MNI space
cit2mni_t1=${out_dir}/CIT168to${out_space}_T1w_head.nii.gz
cit2mni_t2=${out_dir}/CIT168to${out_space}_T2w_head.nii.gz

# ANTs options
nthreads=8
ants_opts="-d 3 -n ${nthreads} -t b -o ${ants_prefix}"

# Calculate CIT168 to MNI152 SyN warp
if [ ! -s ${cit2mni_warp} ]; then
    echo "Calculating SyN warp from CIT168 to ${out_space}"

    case ${atlas_space} in
        "FSL")
            echo "  FSL : using T1w templates only for warp estimation"
            antsRegistrationSyN.sh ${ants_opts} -f ${mni_t1} -m ${cit_t1} 2>&1 > ${logfile}
            ;;
        "2009c")
            echo "  2009c : using both T1w and T2w templates for warp estimation"
            antsRegistrationSyN.sh ${ants_opts} -f ${mni_t1} -f ${mni_t2} -m ${cit_t1} -m ${cit_t2} 2>&1 > ${logfile}
            ;;
        *)
            echo "* Unknown atlas space - exiting"
            exit
    esac

fi

# Apply warp to CIT168 T1w template
if [ ! -s ${cit2mni_t1} ]; then
    echo "Warping CIT168 T1 to ${out_space} space"
    WarpImageMultiTransform 3 ${cit_t1} ${cit2mni_t1} ${cit2mni_warp} ${cit2mni_affine} --use-BSpline
fi

# Apply warp to CIT168 T2w template
if [ ! -s ${cit2mni_t2} ]; then
    echo "Warping CIT168 T2 to ${out_space} space"
    WarpImageMultiTransform 3 ${cit_t2} ${cit2mni_t2} ${cit2mni_warp} ${cit2mni_affine} --use-BSpline
fi

# Apply warp to CIT168 probabilistic atlas

