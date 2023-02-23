#!/usr/bin/env python3
"""
Create CSV of node coordinates from Harvard-Oxford atlas regions in 1mm MNI152 (FSL) space

AUTHOR
----
Mike Tyszka, Ph.D.

LICENSE
----

MIT License

Copyright (c) 2019 Mike Tyszka

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import os
import numpy as np
import xml.etree.ElementTree as ET
import argparse


def main():

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Create node and edge files from Harvard-Oxford atlasz')
    parser.add_argument('--sparse', default=False, action='store_true', help='Decimate nodes')
    parser.add_argument('--intra', default=False, action='store_true', help='Intrahemispheric edges only')
    parser.add_argument('--homo', default=False, action='store_true', help='Homotopic edges only')

    # Parse command line arguments
    args = parser.parse_args()

    # Construct options string for filenames
    options_str = ''
    if args.sparse:
        options_str += '_sparse'
    if args.intra:
        options_str += '_intra'
    if args.homo:
        options_str += '_homo'

    # Output CSV filename
    fstub = 'HarvardOxford'
    node_fname = fstub + options_str + '.node'
    edge_fname = fstub + options_str + '.edge'

    # Excluded labels
    excludes = ['Left Cerebral White Matter',
                'Left Cerebral Cortex ',
                'Left Lateral Ventricle',
                'Right Cerebral White Matter',
                'Right Cerebral Cortex ',
                'Right Lateral Ventricle',
                'Brain-Stem']

    cort_fname = os.path.join(os.environ['FSLDIR'], 'data', 'atlases', 'HarvardOxford-Cortical.xml')
    cort_xml = ET.parse(cort_fname).getroot()

    rh_cort_nodes = []

    for tag in cort_xml.findall('data/label'):

        name = 'Right ' + tag.text

        index = int(tag.get('index'))

        # Map voxel coordinates to MNI152 space
        x, y, z = vox2mni(tag)

        # Force RH (positive x)
        x = np.abs(x)

        # Enforce minimum distance from midline (3.0 mm)
        x = np.clip(x, 3.0, None)

        print('{:80s} {:2d} ({:0.1f}, {:0.1f}, {:0.1f})'.format(name, index, x, y, z))

        rh_cort_nodes.append([name, index, x, y, z])

    # Synthesize LH seed coordinates by reflection about midline

    lh_cort_nodes = []

    for node in rh_cort_nodes:

        name = node[0].replace('Right', 'Left')
        index = node[1] + len(rh_cort_nodes)

        # Mirror about midline
        x = -node[2]
        y = node[3]
        z = node[4]

        print('{:80s} {:2d} ({:0.1f}, {:0.1f}, {:0.1f})'.format(name, index, x, y, z))

        lh_cort_nodes.append([name, index, x, y, z])

    subcort_fname = os.path.join(os.environ['FSLDIR'], 'data', 'atlases', 'HarvardOxford-Subcortical.xml')
    subcort_xml = ET.parse(subcort_fname).getroot()

    subcort_nodes = []

    for tag in subcort_xml.findall('data/label'):

        name = tag.text

        if name in excludes:

            print('  - {}'.format(name))

        else:

            index = int(tag.get('index')) + 2 * len(rh_cort_nodes)

            # Map voxel coordinates to MNI152 space
            x, y, z = vox2mni(tag)

            print('{:80s} {:2d} ({:0.1f}, {:0.1f}, {:0.1f})'.format(name, index, x, y, z))

            subcort_nodes.append([name, index, x, y, z])

    # Optional node culling
    if args.sparse:
        rh_cort_nodes = rh_cort_nodes[0::5]
        lh_cort_nodes = lh_cort_nodes[0::5]
        subcort_nodes = []

    # Concatenate RH cort, LH cort and subcort nodes
    nodes = rh_cort_nodes + lh_cort_nodes + subcort_nodes

    # Write nodes to CSV file
    print('')
    print('Writing nodes to {}'.format(node_fname))

    with open(node_fname, 'w') as fd:

        for node in nodes:

            node, idx, x, y, z = node

            # Replace whitespace with '.'
            node = node.replace(' ', '.')

            fd.write('{:10.3f}{:10.3f}{:10.3f}{:10d}{:10d} {:s}\n'.format(x, y, z, 1, 1, node))


    # Generate demo random edge weights
    n_nodes = len(nodes)
    n_hemi = int(n_nodes / 2)
    w_edge = np.random.random([n_nodes, n_nodes])

    # Intrahemispheric edges only
    if args.intra:
        w_edge[0:n_hemi, n_hemi:] = 0.0
        w_edge[n_hemi:, 0:n_hemi] = 0.0

    # Homotopic edges only
    if args.homo:
        w_intra = w_edge[0:n_hemi, n_hemi:]
        w_homo = np.diag(w_intra)
        w_edge = np.diag(w_homo, n_hemi)

    # Save edge weight array as CSV
    print('Writing random edges to {}'.format(edge_fname))
    np.savetxt(edge_fname,
               w_edge,
               fmt='%0.3f',
               delimiter=',')


def vox2mni(tag):
    """
    fslorient -getsform HarvardOxford-cort-maxprob-thr0-2mm.nii.gz

    -2   0   0   90
     0   2   0 -126
     0   0   2  -72
     0   0   0    1

    Parameters
    ----------
    tag

    Returns
    -------

    """

    # Empirical adjustments for SurfIce mni152_2009 pial mesh
    sf = 1.05
    z_adjust = 3

    xv = int(tag.get('x'))
    yv = int(tag.get('y'))
    zv = int(tag.get('z'))

    x = sf * (-2.0 * xv + 90)
    y = sf * (2.0 * yv - 126)
    z = sf * (2.0 * zv - 72 + z_adjust)

    return x, y, z

if '__main__' in __name__:
    main()
