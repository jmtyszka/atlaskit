#!/bin/bash
# Create a CIFTI dlabel file from a 4D probabilistic atlas in MNI152 space
#
# DEPENDENCIES:
# - HCP workbench installation
# - HCP Pipeline installation with HCPPIPEDIR set
# - FSL installation with FSLDIR set
#
# AUTHORS : Julien Dubois and Mike Tyszka
# PLACE   : Caltech
# DATES   : 2017-02-01 From Julien's solution

#
# Setup
#

# Check for HCP Pipeline env
if [ "x${HCPPIPEDIR}" == "x" ]; then
  echo "HCPPIPEDIR not set - please install and configure HCP Pipelines"
  exit
fi

# Temporary intermediate filenames
tmp_caret=`mktemp`.nii.gz
tmp_prob_atlas=`mktemp`.nii.gz
tmp_det_atlas=`mktemp`.nii.gz

# FLIRT resampling options
FLIRT_opts="-init ${FSLDIR}/etc/flirtsch/ident.mat -applyxfm -interp sinc -sincwindow hanning"

# Use HCPPIPEDIR to locate Atlas_ROIs.2 mesh template
atlas_mesh_hcp=${HCPPIPEDIR}/global/templates/standard_mesh_atlases/Atlas_ROIs.2.nii.gz

#
# Argument parsing
#

if [ $# -lt 2 ]; then
  echo "USAGE : prob2cifti.sh <4D probabilistic atlas> <CIFTI label definitions>"
  exit 
fi

prob_atlas=$1
label_txt=$2

# Create dlabel filename from prob atlas filename
dlabel_nii=${prob_atlas%%.nii.gz}.dlabel.nii

#
# Main actions
#

# Resample probabilistic atlas to same space as HCP atlas mesh template
echo "Resampling MNI152 probabilistic atlas to HCP mesh space"
flirt -in ${prob_atlas} -out ${tmp_prob_atlas} -ref ${atlas_mesh_hcp} ${FLIRT_opts}

# Convert probabilistic atlas to deterministic by thresholding (p >= 0.25)
echo "Converting probabilistic labels to deterministic labels (p >= 0.25)"
prob2det.sh ${tmp_prob_atlas} ${tmp_det_atlas} 0.25

# Convert deterministic atlas to Caret Nifti extension format with label info
echo "Converting deterministic labels to Caret Nifti extension format"
wb_command -volume-label-import ${tmp_det_atlas} ${label_txt} ${tmp_caret}

# Create CIFTI dlabel file from intermediate atlas
echo "Creating CIFTI dlabel file: ${dlabel_nii}"
wb_command -cifti-create-label ${dlabel_nii} -volume ${tmp_caret} ${atlas_mesh_hcp}

# Clean up
echo "Cleaning up"
rm -rf ${tmp_caret} ${tmp_prob_atlas} ${tmp_det_atlas}
