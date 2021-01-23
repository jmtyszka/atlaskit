#!/bin/bash
# Collect necessary images, surfaces and editing files for one subject to
# pass to a FS editor
#
# USAGE : fs_mkeditpack <FS subject recon dir> <output directory>
#
# AUTHOR : Mike Tyszka
# PLACE  : Caltech
# DATES  : 2021-01-21 JMT From scratch

if [ $# -lt 2 ]
then
    echo "USAGE : fs_mkeditpack <FS subject recon dir> <output directory>"
    exit
fi

subj_dir=$1
subj_basedir=$(basename ${subj_dir})

out_dir=$2

pack_dir=${out_dir}/${subj_basedir}
mkdir -p ${pack_dir}/mri
mkdir -p ${pack_dir}/surf
mkdir -p ${pack_dir}/tmp
mkdir -p ${pack_dir}/outbox

# Volume images
cp ${subj_dir}/mri/T1.mgz ${pack_dir}/mri
cp ${subj_dir}/mri/brainmask.mgz ${pack_dir}/mri
cp ${subj_dir}/mri/wm.mgz ${pack_dir}/mri
cp ${subj_dir}/mri/brain.finalsurfs.mgz ${pack_dir}/mri
cp ${subj_dir}/mri/aseg.mgz ${pack_dir}/mri

# Deface T1
mgz_fname=${pack_dir}/mri/T1.mgz
nii_fname=${mgz_fname/.mgz/.nii.gz}
def_fname=${mgz_fname/.mgz/_defaced.nii.gz}
mri_convert ${mgz_fname} ${nii_fname} > /dev/null
pydeface.py -i ${nii_fname} --overwrite
mri_convert ${def_fname} ${mgz_fname} > /dev/null

# Cleanup after defacing
rm -rf ${pack_dir}/mri/*.nii.gz

# Manually edited brain.finalsurfs if present
bfs=${subj_dir}/mri/brain.finalsurfs.manedit.mgz
if [ -s ${bfs} ]
then
    cp ${bfs} ${pack_dir}/mri
fi

# Surfaces
cp ${subj_dir}/surf/lh.pial ${pack_dir}/surf
cp ${subj_dir}/surf/rh.pial ${pack_dir}/surf
cp ${subj_dir}/surf/lh.white ${pack_dir}/surf
cp ${subj_dir}/surf/rh.white ${pack_dir}/surf

# Control points if present
cps=${subj_dir}/tmp/control.dat
if [ -s ${cps} ]
then
    cp ${cps} ${pack_dir}/tmp
fi
