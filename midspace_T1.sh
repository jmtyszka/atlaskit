#!/usr/bin/env bash
# Build an optimal univariate (T1w only) brain midspace using ANTS.
# This will take from several hours to several days depending on the number of
# individual images, image dimensions and computing resources.
#
# SETUP:
#
# 1. Place all T1w .nii.gz images in a single directory.
# 2. Copy this script to the image directory, make sure it's executable and run it in place.
#
# OUTPUT:
#
# The T1w midspace template will be named MIDSPACE_template0.nii.gz.
#
# AUTHOR : Mike Tyszka, Ph.D.
# PLACE  : Caltech
# DATES  : 2016-12-13 JMT Clone from CIT168 atlas tree
#
# MIT License
#
# Copyright (c) 2016 Mike Tyszka
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

# File list for template (both T1w and T2w)
files=`ls *.nii.gz`

# antsMultivariateTemplateConstruction options (not all used)
# -d 3 (3D volumes)
# -o MIDSPACE_ (output prefix)
# -k 1 (1 modality)
# -c 1 (SGE job manager)
# -n 0 (No bias correction - HCP already corrected)
# -w 1x1 (Modality weights)
# -i 1 (Single iteration)

# Full, four iteration templates without CIT168 seed templates
antsMultivariateTemplateConstruction.sh -d 3 -o MIDSPACE_ -k 1 -c 1 -n 0 ${files}

