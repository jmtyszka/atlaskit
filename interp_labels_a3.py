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
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator
from skimage import measure
import time
from scipy.ndimage.filters import gaussian_filter
import random


def ReduceSlices2Contours(Lsub, slices):
    """
    Detects contours of segmentation in a slices
    
    @param Lsub: Label Volume
    @type Lsub: 3D-np.array of e.g. NIFTI image
    @param slices: coordinates of slices in each axis direction
    @type slices: 3-Tupel
    @return: Volume with Labels reduced to contours in same format as input, list of all detected contorus w/ x,y of points in each contour
    @rtype: list
    """
    new_Lsub = np.zeros_like(Lsub)
    all_contours = []
    for axis in range(3):
        for i in slices[axis][0]:
            if axis == 0:
                myslice = Lsub[i,:,:]
            elif axis == 1:
                myslice = Lsub[:,i,:]
            elif axis == 2:
                myslice = Lsub[:,:,i]
            contours = measure.find_contours(myslice, 0.5,fully_connected='high')
            myslice = np.zeros_like(myslice)
            for contour in contours:
                all_contours.append(contour)
                for p in xrange(len(contour)):
                    x,y = contour[p]
                    myslice[int(x), int(y)] = 1
                if axis == 0:
                    new_Lsub[i,:,:] += myslice
                elif axis == 1:
                    new_Lsub[:,i,:] += myslice
                elif axis == 2:
                    new_Lsub[:,:,i] += myslice
    new_Lsub = (new_Lsub > 0.5).astype(int)
    return(new_Lsub, all_contours)


def EvalSliceDistance(slices):
    """
    Get the median distance between slices
    
    @param slices: coordinates of slices in each axis direction
    @type slices: 3-Tupel
    @return: median slice distance
    @rtype: float
    """
    distances = []
    for a in xrange(len(slices)):
        s = slices[a][0]
        if len(s) > 1:
            dl = list(s[1:len(s)] - s[0:(len(s) - 1)])
            for d in dl:
                distances.append(d)
        else:
            distances.append(3)
    return(np.median(distances))


def MakeSamplePoints(vol, slices, dist):
    """
    add regularly spaced points on slice
    
    @param vol: 3D volume of segmentation
    @type vol: 3D np.array
    @param slices: coordinates of slices
    @type slices: 3-Tupel
    @param dist: distance among neirest neighbors (i.e. along x/y axis, not diagonally)
    @type dist: float
    @return: 3D volumen of segmentation
    @rtype: 3D np.array
    """
    vol_s = np.zeros_like(vol)
    for axis in range(3):
        for i in slices[axis][0]:
            if axis == 0:
                myslice = vol[i,:,:]
            elif axis == 1:
                myslice = vol[:,i,:]
            elif axis == 2:
                myslice = vol[:,:,i]
            xs = np.arange(0, myslice.shape[0], dist)
            ys = np.arange(0, myslice.shape[1], dist)
            myslice_s = np.zeros_like(myslice)
            for x in xs:
                for y in ys: 
                    myslice_s[x,y] = 1
            myslice_s = (myslice + myslice_s) > 1
            if axis == 0:
                vol_s[i,:,:] = myslice_s
            elif axis == 1:
                vol_s[:,i,:] = myslice_s
            elif axis == 2:
                vol_s[:,:,i] = myslice_s
    return(vol_s)


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
    
    
def FindSlices(label, n_slices):
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
    Dmin_default = 0.1
    if n_slices[0] != -1:    
        Dmin = Dmin_default * 2
        Sx = (np.array([]),1)
        while Sx[0].shape[0] < n_slices[0] and Dmin >= 0.01:
            Sx = np.where(Dx > Dmin)
            Dmin -= Dmin * .1
    else:
        Sx = np.where(Dx > Dmin_default)

    if n_slices[1] != -1:    
        Dmin = Dmin_default * 2
        Sy = (np.array([]),1)
        while Sy[0].shape[0] < n_slices[1] and Dmin >= 0.01:
            Sy = np.where(Dy > Dmin)
            Dmin -= Dmin * .1
    else:
        Sy = np.where(Dy > Dmin_default)

    if n_slices[2] != -1:    
        Dmin = Dmin_default * 2
        Sz = (np.array([]),1)
        while Sz[0].shape[0] < n_slices[2] and Dmin >= 0.01:
            Sz = np.where(Dz > Dmin)
            Dmin -= Dmin * .1
    else:
        Sz = np.where(Dz > Dmin_default)

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
        
        # Extract slice and flatten
        vv = vol[x,:,:].reshape(-1,1)
        
        # Slice coordinate mesh with ij indexing
        xm, ym, zm = np.meshgrid(x, yv, zv, indexing='ij', )
        xm, ym, zm = xm.reshape(-1,1), ym.reshape(-1,1), zm.reshape(-1,1)
        new_nodes = np.hstack([xm, ym, zm])

        # Append the nodes and values
        nodes = _safe_append(nodes, new_nodes)
        vals = _safe_append(vals, vv)
    
    for y in Sy[0]:
        
        # Extract slice and flatten
        vv = vol[:,y,:].reshape(-1,1)
        
        # Slice coordinate mesh with ij indexing
        xm, ym, zm = np.meshgrid(xv, y, zv, indexing='ij', )
        xm, ym, zm = xm.reshape(-1,1), ym.reshape(-1,1), zm.reshape(-1,1)
        new_nodes = np.hstack([xm, ym, zm])

        # Append the nodes and values
        nodes = _safe_append(nodes, new_nodes)
        vals = _safe_append(vals, vv)
        
    for z in Sz[0]:
        
        # Extract slice and flatten
        vv = vol[:,:,z].reshape(-1,1)

        # Slice coordinate mesh with ij indexing
        xm, ym, zm = np.meshgrid(xv, yv, z, indexing='ij', )
        xm, ym, zm = xm.reshape(-1,1), ym.reshape(-1,1), zm.reshape(-1,1)
        new_nodes = np.hstack([xm, ym, zm])

        # Append the nodes and values
        nodes = _safe_append(nodes, new_nodes)
        vals = _safe_append(vals, vv)
        
    # Remove duplicate locations
    rr = _unique_rows(nodes)
    nodes = nodes[rr,:]
    vals = vals[rr]
    
    print('  Using %d unique nodes' % vals.size)

    return nodes, vals
    

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
    


def alpha_shape(points, tri, alpha):
    classification = np.zeros(tri.vertices.shape[0])

    for i in xrange(tri.vertices.shape[0]):
        pa = points[tri.vertices[i,0]]
        pb = points[tri.vertices[i,1]]
        pc = points[tri.vertices[i,2]]
        pd = points[tri.vertices[i,3]]

        # a = |x_1 y_1 z_1 1; x_2 y_2 z_2 1; x_3 y_3 z_3 1; x_4 y_4 z_4 1|
        a = np.linalg.det(np.array([
            [pa[0], pa[1], pa[2], 1],
            [pb[0], pb[1], pb[2], 1],
            [pc[0], pc[1], pc[2], 1],
            [pd[0], pd[1], pd[2], 1]
        ]))

        # D=[x_1^2+y_1^2+z_1^2 x_1 y_1 z_1 1; x_2^2+y_2^2+z_2^2 x_2 y_2 z_2 1; x_3^2+y_3^2+z_3^2 x_3 y_3 z_3 1; x_4^2+y_4^2+z_4^2 x_4 y_4 z_4 1] 
        D=np.array([[pow(pa[0],2)+pow(pa[1],2)+pow(pa[2],2), pa[0], pa[1], pa[2], 1],
                    [pow(pb[0],2)+pow(pb[1],2)+pow(pb[2],2), pb[0], pb[1], pb[2], 1],
                    [pow(pc[0],2)+pow(pc[1],2)+pow(pc[2],2), pc[0], pc[1], pc[2], 1],
                    [pow(pd[0],2)+pow(pd[1],2)+pow(pd[2],2), pd[0], pd[1], pd[2], 1]
                ])

        D_x = np.linalg.det(D[:,[0,2,3,4]])
        D_y = np.linalg.det(D[:,[0,1,3,4]]) * -1
        D_z = np.linalg.det(D[:,[0,1,2,4]])

        # c=|x_1^2+y_1^2+z_1^2 x_1 y_1 z_1; x_2^2+y_2^2+z_2^2 x_2 y_2 z_2; x_3^2+y_3^2+z_3^2 x_3 y_3 z_3; x_4^2+y_4^2+z_4^2 x_4 y_4 z_4|. 
        c = np.linalg.det(np.array([
            [pow(pa[0],2)+pow(pa[1],2)+pow(pa[2],2), pa[0], pa[1], pa[2]],
            [pow(pb[0],2)+pow(pb[1],2)+pow(pb[2],2), pb[0], pb[1], pb[2]],
            [pow(pc[0],2)+pow(pc[1],2)+pow(pc[2],2), pc[0], pc[1], pc[2]],
            [pow(pd[0],2)+pow(pd[1],2)+pow(pd[2],2), pd[0], pd[1], pd[2]]
        ]))

        # circumsphere        
        circum_r = np.sqrt(pow(D_x,2) + pow(D_y,2) + pow(D_z,2) - 4 * a * c) / (2 * abs(a))
 
        # Here's the radius filter.
        #print circum_r
        if circum_r < 1.0/alpha:
            classification[i] = 1.0

    return(classification)


def save_to_nifti(vol, bb, hdr_nii, out_fname):
    """
    Save vol to nifti-file
    
    @param vol: volume to save 
    @type vol: 3d np.array
    @param hdr_nii: nifti object with relevant header info
    @type hdr_nii: nibabel nifti object
    @param out_fname: output filename
    @type out_fname: str
    @return: ''
    @rtype: None
    """
    # Save interpolated label volume
    tmp_vol = np.zeros_like(hdr_nii.get_data())
    tmp_vol = InsertSubVol(tmp_vol, vol, bb)
    out_nii = nib.Nifti1Image(tmp_vol, hdr_nii.get_affine())
    out_nii.to_filename(out_fname)


def SetValsPoints(points, vals, hdr_nii):
    """
    Construct a 3D array with vals at points
    
    @param points: coordinates at which to write values
    @type points: list of tupels
    @param vals: values to write at coordinates
    @type vals: list
    @param hdr_nii: nifti file with relevant header info
    @type hdr_nii: nibabel nifti file object
    @return: 3D np.array with vals at points
    @rtype: np.array
    """
    tmp_vol = np.zeros_like(hdr_nii.get_data())
    #print new_labels.shape
    for i in xrange(len(vals)):
        point = points[i,:]
        tmp_vol[point[0], point[1], point[2]] = vals[i]

    return tmp_vol


def smooth_labels(vol):
    """
    Smooth labels with Gaussian filter
    
    @param vol: label image
    @type vol: 3d np.array
    @return: smoothed label image
    @rtype: 3d np.array
    """
    # Smooth target label region
    print('  Gaussian smoothing original target label')
    vol = gaussian_filter(vol.astype(float), sigma=1.0)
    
    # Normalize smoothed intensities
    vol = vol / vol.max()

    # Threshold smoothed mask at 0.5 to create new boolean mask
    print('  Thresholding smoothed label')
    vol = vol > 0.5


def main():

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Interpolate labels')
    parser.add_argument('-i','--input', required=True, help="Labeled volume")
    parser.add_argument('-l','--labels', help="Label numbers to interpolate, separated by comma")
    parser.add_argument('-p', '--save-preproc', help="Save result of preprocessing", default=False, action='store_const', const=True, dest='save_preproc')
    parser.add_argument('-d', '--save-delaunay', help="Save result of Delaunay tesselation", default=False, action='store_const', const=True, dest='save_delaunay')
    parser.add_argument('-s', '--smooth-results', help="Smooth results of interpolation", default=False, action='store_const', const=True, dest='smooth_labels')
    parser.add_argument('-sl','--slices', help="Label numbers to interpolate, separated by comma")

    # Parse command line arguments
    args = parser.parse_args()

    # Get mandatory filename argument
    label_fname = args.input
    print label_fname

    if args.labels:
        sink = args.labels
        sink = sink.split(',')
        label_nos = []
        for i in xrange(len(sink)):
            label_nos.append(int(sink[i]))
    else:
        # Construct list of unique label values in image
        label_nos = np.unique(labels)

    if args.slices:
        sink = args.slices
        sink = sink.split(',')
        n_slices = []
        for i in xrange(len(sink)):
            n_slices.append(int(sink[i]))
    else:
        # Construct list of unique label values in image
        n_slices = [-1,-1,-1]

    # loop over each unique label value
    for label in label_nos:
        
        if label > 0:
            
            print('Interpolating label %d' % label)

    # Construct output filename    
    out_stub, out_ext = os.path.splitext(label_fname)
    if out_ext == '.gz':
        out_stub, _ = os.path.splitext(out_stub)
    
    # Load labeled volume
    label_nii = nib.load(label_fname)
    labels = label_nii.get_data()

    # Size of image space
    nx, ny, nz = labels.shape

    # Remember when we started processing input
    start_time = time.time()

    # Extract current label
    L = (labels == label).astype(float)

    # Extract minimum subvolume containing label
    Lsub, bb = ExtractMinVol(L)

    # Detect slices in segmentaion image
    slices = FindSlices(Lsub, n_slices)
    print "Number of slices, x: %s, y: %s, z: %s" % (slices[0][0].shape[0],slices[1][0].shape[0],slices[2][0].shape[0])

    # Reduce to segmentation within slices to contour lines to speed up processing
    Lsub_contour, contour = ReduceSlices2Contours(Lsub, slices)

    # Calc median distance of slices
    dist = EvalSliceDistance(slices)

    # Set alpha so that it tetrahedrons between slices are not removed by alpha shape approach
    alpha = 1.0/dist/2

    # Fill contours with sample points for Delaunay tesselation
    Lsub_sub = MakeSamplePoints(Lsub, slices, dist)

    # Input to alpha shapes is the combination of contours and sample points
    Lsub = Lsub_contour + Lsub_sub
    Lsub = (Lsub > 0.5).astype(int)

    if args.save_preproc:
        print('Saving result of preprocessing to %s' % (out_stub + '_preproc.nii.gz'))
        save_to_nifti(Lsub, bb, label_nii, out_stub + '_preproc.nii.gz')

    print('  Label contains %d voxels' % np.sum(Lsub[:]))

    # get coordinates of segmentation labels
    x,y,z = np.where(Lsub > 0)
    points = np.transpose(np.array((x,y,z)))

    # perform Delaunay tesselation
    tri = Delaunay(points)

    # Construct interpolation mesh for volume
    nx, ny, nz = Lsub.shape
    xv, yv, zv = np.arange(0,nx), np.arange(0,ny), np.arange(0,nz)
    xi, yi, zi = np.meshgrid(xv, yv, zv, indexing='ij')
    xi, yi, zi = xi.reshape(-1,1), yi.reshape(-1,1), zi.reshape(-1,1)
    
    new_points = np.zeros((xi.shape[0], 3))
    for i in xrange(xi.shape[0]):
        new_points[i,:] = xi[i], yi[i], zi[i]

    # determine for each point in which tetrahedron it is
    simplices_i = tri.find_simplex(new_points)

    if args.save_delaunay:
        vals = simplices_i.copy()
        vals[vals == -1] = 0.0
        tmp_vol = SetValsPoints(new_points, vals, label_nii)
        print('Saving result of Delaunay tesselation to %s' % (out_stub + '_delaunay.nii.gz'))
        save_to_nifti(tmp_vol, bb, label_nii, out_stub + '_delaunay.nii.gz')


    # perform alpha shape 3 
    print "Vertices in Delaunay tesselation: %s" % tri.vertices.shape[0]
    v_class = alpha_shape(points, tri, alpha)
    print "Vertices in Alpha Complex: %s" % np.sum(v_class)


    print "Points contained in Dalauny tesselation: %s" % len(np.where(simplices_i > -1)[0])
    for i in xrange(len(simplices_i)):
        if simplices_i[i] > -1:
            if v_class[simplices_i[i]] == 0:
                simplices_i[i] = -1
    print "Points contained in alpha complex: %s" % len(np.where(simplices_i > -1)[0])

    # create segmentation image of interpolation
    vals = np.zeros_like(simplices_i)
    vals[simplices_i > -1] = label
    voli = SetValsPoints(new_points, vals, label_nii)

    # smooth labels
    if args.smooth_labels:
        smooth_labels(voli)

    print('Saving result of interpolation to %s' % (out_stub + '_interp.nii.gz'))
    save_to_nifti(voli, bb, label_nii, out_stub + '_interp.nii.gz')

    print('Total Processing Time: %s s' % (time.time() - start_time))

    # Clean exit
    sys.exit(0)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
