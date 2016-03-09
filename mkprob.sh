#!/bin/bash
# Create probabilistic atlases of all subdivisions and the whole amygdala
#
# Labels 1..10 inclusive are amygdala subregions
#
# AUTHOR : Mike Tyszka
# PLACE  : Caltech
# DATES  : 2015-07-30 JMT From scratch

root_dir=..
atlas_dir=${root_dir}/atlas
val_dir=${root_dir}/validation

# Key filenames
pAll_LH=${atlas_dir}/pAll_LH.nii.gz
pAmy_LH=${atlas_dir}/pAmy_LH.nii.gz
pAmyOnly_LH=${atlas_dir}/pAmyOnly_LH.nii.gz

echo ""
echo "Constructing amygdala probabilistic atlas"
echo "----"
probabilistic.py -o ${atlas_dir}/pAll_LH.nii.gz ${val_dir}/*/*_all_*.nii.gz

# Crop atlas to first ten labels
echo ""
echo "Isolating amygdala subdivisions labels"
echo "----"
fslroi ${pAll_LH} ${pAmy_LH} 0 -1 0 -1 0 -1 0 10

# Collapse 4th dimension to generate whole amygdala prob atlas
echo ""
echo "Collapsing all labels into whole amygdala label"
echo "----"
fslmaths ${pAmy_LH} -Tmean -mul 10.0 ${pAmyOnly_LH}

# Add whole amy as final volume in prob atlas
# Atlas now includes all individual subdivisions and a final, whole amygdala label
echo
echo "Adding whole amygdala label to atlas"
echo "----"
fslmerge -t ${pAmy_LH} ${pAmy_LH} ${pAmyOnly_LH}

