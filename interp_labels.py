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
>>> interp_labels.py labels.nii.gz 1 3 4

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
from scipy.spatial import ConvexHull, Delaunay
from scipy.signal import medfilt

def main():
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Interpolate labels')
    parser.add_argument('-i','--input', help="Labeled volume")
    parser.add_argument('-l','--labels', help="Label numbers to interpolate")

    # Parse command line arguments
    args = parser.parse_args()
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
    new_labels = np.zeros_like(labels, dtype='float')
    
    if args.labels:
        label_nos = args.labels
    else:
        # Construct list of unique label values in image
        label_nos = np.unique(labels)
        
    # loop over each unique label value
    for label in label_nos:
        
        if label > 0:
            
            print('Interpolating label %d' % label)

            # Extract current label
            L = (labels == label)
            
            # Extract minimum subvolume containing label
            Lsub, bb = ExtractMinVol(L)
            
            # Find locations of single labeled slices in each axis
            slices = FindSlices(Lsub)
            
            # Construct point value lists over all slices
            nodes, vals = InsideOutside(Lsub, slices)
            
            # RBF Interpolate values within subvolume
            Lsubi = RBFInterpolate(Lsub, nodes, vals)

    
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
    
    bb = (xmin, xmax), (ymin, ymax), (zmin, zmax)
    
    subvol = label[xmin:xmax, ymin:ymax, zmin:zmax]
    
    return subvol, bb
    
    
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
    
    
def InsideOutside3D(vol, slices):
    '''
    '''
    
    nodes = np.array([])
    vals = np.array([])

    # X, Y and Z slice lists
    Sx, Sy, Sz = slices
    
    # Volume dimensions
    nx, ny, nz = vol.shape
    xv, yv, zv = np.arange(0,nx), np.arange(0,ny), np.arange(0,nz)
    
    for x in Sx[0]:
        
        # Extract slice
        Ixy = vol[x,:,:]
    
    for y in Sy[0]:
        
        # Extract slice
        Ixy = vol[:,y,:]

    for z in Sz[0]:
        
        # Extract slice
        Ixy = vol[:,:,z]

        # Slice coordinate mesh        
        xym = np.meshgrid(xv, yv, indexing='ij')
        nodes_xy = np.array(xym).transpose().reshape(-1,3)



    # Append the nodes and values
    nodes = np.append(nodes, nodes_xy, axis=0)
    vals = np.append(vals, vals_xy, axis=0)
        

    return nodes, vals
    

def InsideOutside2D(img):
    '''
    '''
    
    # Generate whole slice coordinate mesh


def RBFInterpolate(vol, nodes, vals):
    '''
    Interpolate binary values at nodes over the volume provided
    '''
    
    voli = np.zeros_like(vol)
    
    return voli
    

# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
    
    
#
# Archive of older ideas
#

#            # Center of mass of hull
#            hx = x[hull.vertices]
#            hy = y[hull.vertices]
#            hz = z[hull.vertices]
#            cx, cy, cz = hx.mean(), hy.mean(), hz.mean()
#            
#            # Inflation scale factor for each node
#            dx, dy, dz = hx-cx, hy-cy, hz-cz
#            r = np.sqrt(dx**2 + dy**2 + dz**2)
#            sf = (r + dr) / r
#            
#            # Inflate hull
#            ihx = cx + dx * sf
#            ihy = cy + dy * sf
#            ihz = cz + dz * sf
#            
#            # Ensure all coordinates are within the image space
#            ihx = np.clip(ihx, 0, nx-1)
#            ihy = np.clip(ihy, 0, ny-1)
#            ihz = np.clip(ihz, 0, nz-1)
#            
#            # Boundary value on inflated hull
#            ihd = np.zeros_like(ihx, dtype='float')
#            
#            # Create N x 3 array of inflated hull points
#            ihp = np.array([ihx,ihy,ihz]).transpose().reshape(-1,3)
#            
#            # Append boundary conditions
#            x = np.append(x, ihx)
#            y = np.append(y, ihy)
#            z = np.append(z, ihz)
#            d = np.append(d, ihd)
#    
#            # Construct RBF interpolator
#            print('  Constructing RBF interpolant')
#            rbf = Rbf(x, y, z, d, function='multiquadric')
#            
#            # Construct interpolation mesh within inflated hull
#            xmin, xmax = np.ceil(ihx.min()), np.ceil(ihx.max())
#            ymin, ymax = np.ceil(ihy.min()), np.ceil(ihy.max())
#            zmin, zmax = np.ceil(ihz.min()), np.ceil(ihz.max())
#            xiv = np.arange(xmin, xmax)
#            yiv = np.arange(ymin, ymax)
#            ziv = np.arange(zmin, zmax)
#            pim = np.meshgrid(xiv, yiv, ziv, indexing='ij')
#            pim = np.array(pim).transpose().reshape(-1,3)
#            
#            # Create mask for points within inflated hull
#            del_hull = Delaunay(ihp)
#            mask = del_hull.find_simplex(pim) >= 0
#            
#            # Extract x, y, z coordinates of mesh points inside inflated hull
#            xi, yi, zi = pim[mask,0], pim[mask,1], pim[mask,2]
#            
#            print('  Interpolating label')
#            Li = rbf(xi,yi,zi)
#            
#            # Update label volume
#            new_labels[xi.astype(int), yi.astype(int), zi.astype(int)] = Li