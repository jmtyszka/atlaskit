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
from scipy.spatial import ConvexHull

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
            mask = (labels == label)
            
            # Find coordinates of all True voxels (n x 3)
            print('  Locating label voxels')
            points = np.where(mask)
            x, y, z = points
            
            # Create bounding box for label
            # Make upper bound inclusive
            xmin, xmax = x.min(), x.max()+1
            ymin, ymax = y.min(), y.max()+1
            zmin, zmax = z.min(), z.max()+1
            print('  Bounding box : (%d, %d, %d) : (%d, %d, %d)' % (xmin, ymin, zmin, xmax, ymax, zmax))
            
            # Construct voxel coordinate mesh within BB. Note ij indexing to match np.where, etc
            print('  Creating interpolation mesh')
            xx, yy, zz = np.arange(xmin, xmax), np.arange(ymin, ymax), np.arange(zmin, zmax)
            xi,yi,zi = np.meshgrid(xx,yy,zz, indexing='ij')
            
            # Count number of voxels in label
            nvox = np.sum(mask)
            print('  Label volume : %d voxels' % nvox)
            
            # All node values are one
            d = np.ones_like(x, dtype='float')
            
            # Find convex hull
            hull = ConvexHull(points)
            
            # Center of mass of hull
            hx = points[hull.vertices,0]
            hy = points[hull.vertices,1]
            hz = points[hull.vertices,2]
            cx, cy, cz = hx.mean(), hy.mean(), hz.mean()
            dx, dy, dz = hx-cx, hy-cy, hz-cz
            
            # Inflate hull by 10%
            hxi = hx + dx * 1.1
            hyi = hy + dy * 1.1
            hzi = hz + dz * 1.1
            
            # Append boundary conditions
            x = np.append(x, hxi)
            y = np.append(y, hyi)
            z = np.append(z, hzi)
            d = np.append(d, np.zeros_like(hxi))
            
            # RBF interpolate
            print('  Constructing RBF interpolant')
            rbf = Rbf(x, y, z, d)
    
            # Interpolate label only within BB
            print('  Interpolating label')
            Li = rbf(xi,yi,zi)
            
            # Update label volume
            new_labels[xmin:xmax,ymin:ymax,zmin:zmax] = Li
            
    
    # Save interpolated label volume
    print('Saving interpolated labels to %s' % out_fname)
    out_nii = nib.Nifti1Image(new_labels, label_nii.get_affine())
    out_nii.to_filename(out_fname)
        
    
    # Clean exit
    sys.exit(0)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()