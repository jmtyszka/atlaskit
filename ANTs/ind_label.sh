#!/bin/bash
#
# Warp atlas labels from master (350um) template to individual spaces (700um)
# Note : this is the script for the master atlas group (F26-30)
# 
# AUTHOR : Mike Tyszka
# PLACE  : Caltech
# DATES  : 2015-03-25 JMT From scratch
#          2015-04-14 JMT Switch to label warping
#          2015-05-06 JMT Correct transform labeling (I2M, etc)
#          2015-05-12 JMT Move to anti-aliased ML interpolation
#
# Copyright 2015 California Institute of Technology
# All rights reserved

# Root directory for all templates/atlases
ROOT_DIR=/Users/jmt/HCP/Atlases

# Master atlas directories and files
MASTER_DIR=${ROOT_DIR}/F26-30_DoubleT1T2
MASTER_ATLAS_DIR=${MASTER_DIR}/ATLAS/current
MASTER_MIDSPACE_DIR=${MASTER_DIR}/MIDSPACE
MASTER_MORPH_DIR=${MASTER_DIR}/MORPHOMETRY
MASTER_ATLAS=${MASTER_ATLAS_DIR}/Caltech_F26-30_Atlas_Labels_350um_BH.nii.gz

# Nearest neighbor interpolation for atlas labels
INTERP=--use-NN

if [ ! -d ${MASTER_MORPH_DIR} ]; then
  mkdir -p ${MASTER_MORPH_DIR}
fi

# Loop over all T1w brains in master midspace folder
for IND_BRAIN in ${MASTER_MIDSPACE_DIR}/??????_T1w_brain_trim.nii.gz
do

  # Extract 6-digit HCP subject ID
  tmp=`basename $IND_BRAIN`
  SID=${tmp%%_T1w_brain_trim.nii.gz}

  echo "Processing subject ${SID}"

  # Key transform files
  # M : Master template space
  # I : Individual space
  I2M_AFFINE=`ls ${MASTER_MIDSPACE_DIR}/M${SID}*Affine.txt`
  I2M_WARP=`ls ${MASTER_MIDSPACE_DIR}/M${SID}*Warp.nii.gz | tail -n 1`
  M2I_WARP=`ls ${MASTER_MIDSPACE_DIR}/M${SID}*InverseWarp.nii.gz`

  # Trimmed T1w brain (for brain mask)
  T1=${MASTER_MIDSPACE_DIR}/${SID}_T1w_brain_trim.nii.gz

  # Individual atlas and brain mask filenames
  IND_ATLAS=${MASTER_MORPH_DIR}/${SID}_Atlas.nii.gz
  IND_BRAIN_MASK=${MASTER_MORPH_DIR}/${SID}_BrainMask.nii.gz

  # Warp atlas to individual space with nearest neighbor interpolation
  if [ ! -s ${IND_ATLAS} ]; then
    echo "  Warping atlas to individual space"
    # Note order of transforms: inverse affine, then inverse warp
    M2I_TRANSFORM="-i ${I2M_AFFINE} ${M2I_WARP}"
    WarpImageMultiTransform 3 ${MASTER_ATLAS} ${IND_ATLAS} -R ${IND_BRAIN} ${INTERP} ${M2I_TRANSFORM} # > /dev/null 2>&1
  fi

  # Create individual brain mask from trimmed T1w image (simple > 0 binarization)
  if [ ! -s ${IND_BRAIN_MASK} ]; then
    echo "  Creating brain mask"
    fslmaths ${T1} -bin ${IND_BRAIN_MASK}
  fi

done

echo "Done - all output in ${MASTER_MORPH_DIR}"
