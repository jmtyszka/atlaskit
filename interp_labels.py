#!/opt/local/bin/python
'''
Interpolate label between sparse sections
- Speeds up manual labeling for larger label volumes

Usage
----
interp_labels.py <label image> <index list>
interp_labels.py -h

Example
----
>>> interp_labels.py -i labels.nii.gz -l 1,3,4

Authors
----
Mike Tyszka, Caltech Brain Imaging Center

Dates
----
2015-09-28 JMT From scratch

License
----
This file is part of atlaskit.

    atlaskit is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    atlaskit is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with atlaskit.  If not, see <http://www.gnu.org/licenses/>.

Copyright
----
2015 California Institute of Technology.
'''

__version__ = '0.1.0'

import os
import sys
import argparse
import nibabel as nib
import numpy as np
from scipy.interpolate import Rbf
from scipy.signal import medfilt
from scipy.ndimage.morphology import distance_transform_edt as DT


def main():
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Interpolate labels')
    parser.add_argument('-i','--input', required=True, help="Labeled volume")
    parser.add_argument('-l','--labels', help="Label numbers to interpolate, separated by comma")

    # Parse command line arguments
    args = parser.parse_args()

    # Get mandatory filename argument
    label_fname = args.input

    # Construct output filename    
    out_stub, out_ext = os.path.splitext(label_fname)
    if out_ext == '.gz':
        out_stub, _ = os.path.splitext(out_stub)
    out_fname = out_stub + '_interp.nii.gz'

    # Load labeled volume
    label_nii = nib.load(label_fname)
    labels = label_nii.get_data()
    
    # Size of image space
    nx, ny, nz = labels.shape
    
    # Destination label volume
    new_labels = labels.copy()
    
    if args.labels:
        sink = args.labels
        sink = sink.split(',')
        label_nos = []
        for i in xrange(len(sink)):
            label_nos.append(int(sink[i]))
    else:
        # Construct list of unique label values in image
        label_nos = np.unique(labels)

    # loop over each unique label value
    for label in label_nos:
        
        if label > 0:
            
            print('Interpolating label %d' % label)

            # Extract current label
            L = (labels == label).astype(float)
            
            # Extract minimum subvolume containing label
            Lsub, bb = ExtractMinVol(L)
            
            print('  Label contains %d voxels' % np.sum(Lsub[:]))

            # Find locations of single labeled slices in each axis
            slices = FindSlices(Lsub)
            
            # Count slices
            nSx = slices[0][0].size
            nSy = slices[1][0].size
            nSz = slices[2][0].size
            
            # Report number of slices detected
            print('  X slices : %d' % nSx)
            print('  Y slices : %d' % nSy)
            print('  Z slices : %d' % nSz)
            
            # Only interpolate if slice-like features found
            if nSx > 1 or nSy > 1 or nSz > 1:
            
                # Construct point value lists over all slices
                nodes, vals = NodeValues(Lsub, slices)
                
                # RBF Interpolate values within subvolume
                # Returns thresholded integer volume
                Lsubi = RBFInterpolate(Lsub, nodes, vals)
                
                # Scale interpolation back to original label value
                Lsubi *= label
                
                # Insert interpolated volume back into new label volume
                new_labels = InsertSubVol(new_labels, Lsubi, bb)
                
            else:
                
                print('Insufficient slice-like features found - skipping label')

    
    # Save interpolated label volume
    print('Saving interpolated labels to %s' % out_fname)
    out_nii = nib.Nifti1Image(new_labels, label_nii.get_affine())
    out_nii.to_filename(out_fname)
        
    
    # Clean exit
    sys.exit(0)

    
def ExtractMinVol(label):
    '''
    Extract minimum subvolume containing label voxels
    '''
    
    ii = np.argwhere(label)    
    
    xmin, xmax = ii[:,0].min(), ii[:,0].max()+1
    ymin, ymax = ii[:,1].min(), ii[:,1].max()+1
    zmin, zmax = ii[:,2].min(), ii[:,2].max()+1
    
    bb = xmin, xmax, ymin, ymax, zmin, zmax
    
    subvol = label[xmin:xmax, ymin:ymax, zmin:zmax]
    
    return subvol, bb

 
def InsertSubVol(label, new_subvol, bb):
    '''
    Insert subvolume at a given location in label volume
    Treat zeros in new subvolume as transparent to avoid unnecessary
    overwriting of other labels in the original label volume
    '''
    
    # Unpack subvolume limits from bounding box
    xmin, xmax, ymin, ymax, zmin, zmax = bb

    # Find subvolume dimensions
    nx, ny, nz = new_subvol.shape
    
    # Extract corresponding subvolume from original label volume
    label_subvol = label[xmin:(xmin+nx), ymin:(ymin+ny), zmin:(zmin+nz)]
    
    # Create alpha mask for new subvolume
    subvol_mask = new_subvol > 0.0
    
    # Combine old and new subvolumes
    label_subvol[subvol_mask] = new_subvol[subvol_mask]

    # Insert modified subvolume into original label volume
    label[xmin:(xmin+nx), ymin:(ymin+ny), zmin:(zmin+nz)] = label_subvol
    
    return label
    
    
def FindSlices(label):
    '''
    Locate likely isolated slices in each axis
    '''
    
    # Integral over volume
    ii = float(label.ravel().sum())
    
    # Project onto each axis
    Px = np.sum(np.sum(label,axis=2), axis=1) / ii
    Py = np.sum(np.sum(label,axis=2), axis=0) / ii
    Pz = np.sum(np.sum(label,axis=1), axis=0) / ii
    
    # Subtract median filtered baseline estimate (k = 3)
    # Retains only spikes
    Dx = Px - medfilt(Px)
    Dy = Py - medfilt(Py)
    Dz = Pz - medfilt(Pz)
    
    # Locate spikes in residual
    Dmin = 0.1
    Sx = np.where(Dx > Dmin)
    Sy = np.where(Dy > Dmin)
    Sz = np.where(Dz > Dmin)
    
    return Sx, Sy, Sz
    
    
def NodeValues(vol, slices):
    '''
    Generate coordinate nodes and values for interpolation
    Extract values from x, y and z slices in volume
    '''
    
    # Init coordinate and value arrays
    nodes = np.array([])
    vals = np.array([])

    # X, Y and Z slice lists
    Sx, Sy, Sz = slices
    
    # Volume dimensions
    nx, ny, nz = vol.shape
    xv, yv, zv = np.arange(0,nx), np.arange(0,ny), np.arange(0,nz)
    
    for x in Sx[0]:
        
        # Extract slice and remove singlet dimension
        v_yz = vol[x,:,:].squeeze()
        
        # Inside-outside values and nodes from slice
        io_yz, yy, zz = InsideOutside(v_yz)
        
        # Append the nodes and values
        nodes = _safe_append(nodes, new_nodes)
        vals = _safe_append(vals, vv)
    
    for y in Sy[0]:
        
        # Do nothing for now
        
    for z in Sz[0]:

        # Do nothing for now        

        
    # Remove duplicate locations
    rr = _unique_rows(nodes)
    nodes = nodes[rr,:]
    vals = vals[rr]
    
    print('  Using %d unique nodes' % vals.size)

    return nodes, vals
    
def InsideOutside(s):
    '''
    Create inside-outside function for slice and extract nodes, values
    just inside, just outside and on the boundary
    
    Arguments
    ----
    s : 2D numpy integer array
        Extracted slice of label volume
    '''
    
    nx, ny = s.shape

    # Boundary voxel mask (Canny edge detection)
    
    # Inside-outside distance transform
    # Positive outside, negative inside
    d = DT(1-s) - DT(s)
    
    return d
    
    

def RBFInterpolate(vol, nodes, vals, function='multiquadric', smooth=0.5):
    '''
    Interpolate node values within the volume using a radial basis function
    '''

    # Construct RBF interpolator from node values
    print('  Constructing interpolator')
    print('    Function  : %s' % function)
    print('    Smoothing : %0.1f' % smooth)
    xx, yy, zz, vv = nodes[:,0], nodes[:,1], nodes[:,2], vals[:]
    rbf = Rbf(xx, yy, zz, vv, function=function, smooth=smooth)   
    
    # Construct interpolation mesh for volume
    nx, ny, nz = vol.shape
    xv, yv, zv = np.arange(0,nx), np.arange(0,ny), np.arange(0,nz)
    xi, yi, zi = np.meshgrid(xv, yv, zv, indexing='ij')
    xi, yi, zi = xi.reshape(-1,1), yi.reshape(-1,1), zi.reshape(-1,1)

    # Interpolate over entire volume
    print('  Interpolating subvolume over %d voxels' % xi.size)
    voli = rbf(xi, yi, zi).reshape(nx,ny,nz)
    
    # Threshold at 0.5
    voli = (voli > 0.5).astype(int)
    
    return voli

    
def _safe_append(aa, bb, axis=0):
    '''
    Append new_nodes to nodes along the given axis
    Allow for nodes = [], in which case just new_nodes is returned
    Remove and duplicate rows before returrning
    '''

    if not aa.any():
        aabb = bb.copy()
    else:
        aabb = np.append(aa, bb, axis=axis)
    
    return aabb

    
def _unique_rows(a):
    '''
    Return logical mask of unique rows in 2D array
    Based on best answer to http://stackoverflow.com/questions/8560440
    '''
    order = np.lexsort(a.T)
    a = a[order]
    diff = np.diff(a, axis=0)
    ui = np.ones(len(a), 'bool')
    ui[1:] = (diff != 0).any(axis=1)
    
    return ui
    

# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
