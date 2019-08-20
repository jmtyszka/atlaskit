#!/bin/bash
# ANTs B-spline SyN warp CIT168 head template(s) to an MNI space
# - currently supports FSL MNI152 1 mm or MNI ICBM152 asym nonlin 2009c 1 mm spaces
# - uses whole head templates for warping
# - T2w templates used if available (eg with 2009c)
# - Gray matter probabilistic tissue mask
#
# USAGES
# AUTHOR : Mike Tyszka
# PLACE  : Caltech
# DATES  : 2017-10-31 JMT Adapt from mirror_warp_calc.sh
#          2018-05-30 JMT Generalize to other MNI spaces
#          2019-08-20 JMT Add gm mask
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
# 2018-2019 California Institute of Technology.

# Key directories (edit as needed)
root_dir="${HOME}/Box/CIT168_SubCortical_Atlas"
template_dir="${root_dir}/Templates"
icmb_dir="${HOME}/Box/Atlases/MNI_ICBM152_2009c/mni_icbm152_nlin_asym_09c"

# CIT168 700 um templates and tissue masks
cit_t1_brain=${template_dir}/CIT168_T1w_700um.nii.gz
cit_t2_brain=${template_dir}/CIT168_T2w_700um.nii.gz
cit_t1_head=${template_dir}/CIT168_T1w_head_700um.nii.gz
cit_t2_head=${template_dir}/CIT168_T2w_head_700um.nii.gz
cit_gm=${template_dir}/CIT168_gm_700um.nii.gz

# Reinforcement learning prob atlas
cit_prob="${root_dir}/Labeling/reinf_learn_atlas/prob_atlas_bilateral.nii.gz"

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

# MNI152 T1w and T2w whole head templates
case ${atlas_space} in
    "FSL")
        out_space="MNI152-FSL"
        mni_t1_head=${FSLDIR}/data/standard/MNI152_T1_1mm.nii.gz
        mni_t2_head=""
        ;;
    "2009c")
        out_space="MNI152-2009c"
        mni_t1_head=${icmb_dir}/mni_icbm152_t1_tal_nlin_asym_09c.nii
        mni_t2_head=${icmb_dir}/mni_icbm152_t2_tal_nlin_asym_09c.nii
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

# ANTs options
nthreads=8
ants_opts="-d 3 -n ${nthreads} -t b -o ${ants_prefix}"

# Calculate CIT168 to MNI152 SyN warp
if [ ! -s ${cit2mni_warp} ]; then
    echo "Calculating SyN warp from CIT168 to ${out_space}"

    case ${atlas_space} in
        "FSL")
            echo "  FSL : using T1w templates only for warp estimation"
            antsRegistrationSyN.sh ${ants_opts} -f ${mni_t1_head} -m ${cit_t1_head} 2>&1 > ${logfile}
            ;;
        "2009c")
            echo "  2009c : using both T1w and T2w templates for warp estimation"
            antsRegistrationSyN.sh ${ants_opts} -f ${mni_t1_head} -f ${mni_t2_head} -m ${cit_t1_head} -m ${cit_t2_head} 2>&1 > ${logfile}
            ;;
        *)
            echo "* Unknown atlas space - exiting"
            exit
    esac

fi

# Apply CIT168 to MNI warp to T1w, T2w (brain and head) and prob atlas images

echo "Warping CIT168 images to ${out_space} space"

# CIT168 templates in MNI space
cit2mni_t1_brain=${out_dir}/CIT168to${out_space}_T1w_brain.nii.gz
cit2mni_t2_brain=${out_dir}/CIT168to${out_space}_T2w_brain.nii.gz
cit2mni_t1_head=${out_dir}/CIT168to${out_space}_T1w.nii.gz
cit2mni_t2_head=${out_dir}/CIT168to${out_space}_T2w.nii.gz
cit2mni_gm=${out_dir}/CIT168to${out_space}_gm.nii.gz
cit2mni_prob=${out_dir}/CIT168to${out_space}_prob.nii.gz
cit2mni_det=${out_dir}/CIT168to${out_space}_det.nii.gz

echo "----------------"
echo "T1w brain"
WarpImageMultiTransform 3 ${cit_t1_brain} ${cit2mni_t1_brain} -R ${mni_t1_head} ${cit2mni_warp} ${cit2mni_affine} --use-BSpline

echo "----------------"
echo "T2w brain"
WarpImageMultiTransform 3 ${cit_t2_brain} ${cit2mni_t2_brain} -R ${mni_t1_head} ${cit2mni_warp} ${cit2mni_affine} --use-BSpline

echo "----------------"
echo "T1w head"
WarpImageMultiTransform 3 ${cit_t1_head} ${cit2mni_t1_head} -R ${mni_t1_head} ${cit2mni_warp} ${cit2mni_affine} --use-BSpline

echo "----------------"
echo "T2w head"
WarpImageMultiTransform 3 ${cit_t2_head} ${cit2mni_t2_head} -R ${mni_t1_head} ${cit2mni_warp} ${cit2mni_affine} --use-BSpline

echo "----------------"
echo "Gray matter mask"
WarpImageMultiTransform 3 ${cit_gm} ${cit2mni_gm} -R ${mni_t1_head} ${cit2mni_warp} ${cit2mni_affine} --use-BSpline

echo "----------------"
echo "Probabilistic labels"
WarpTimeSeriesImageMultiTransform 4 ${cit_prob} ${cit2mni_prob} -R ${mni_t1_head} ${cit2mni_warp} ${cit2mni_affine}

# Calculate warped deterministic labels from warp probabilistic labels (p >= 0.25)
echo "----------------"
echo "Deterministic labels (p > 0.25)"
prob2det.sh ${cit2mni_prob} ${cit2mni_det} 0.25

echo "Done"
