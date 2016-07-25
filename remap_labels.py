#!/usr/bin/env python3
"""
Use a new key to reorder label numbers in label volumes.

Usage
----
remap_labels.py -ok <Old Key> -nk <New Key> <Label Volumes> ...
remap_labels.py -h

Example
----
>>> remap_labels.py -ok old.txt -nk new.txt labels1.nii.gz labels2.nii.gz

Authors
----
Mike Tyszka, Caltech, Division of Humaninities and Social Sciences

Dates
----
2015-11-19 JMT From scratch

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
"""

__version__ = '0.1.0'

import os, sys
import argparse
import nibabel as nib
import numpy as np
import pandas as pd


def main():
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Reorder labels using a new key')
    parser.add_argument('-ok','--oldkey', help="Old ITK-SNAP label key file")
    parser.add_argument('-nk','--newkey', help="New ITK-SNAP label key file")
    parser.add_argument('labels', metavar='N', type=str, nargs='+',
                        help='list of label volume files to reorder')
    
    # Parse command line arguments
    args = parser.parse_args()
    old_key_fname = args.oldkey
    new_key_fname = args.newkey
    label_fnames = args.labels
    
    # Load old label key
    if os.path.isfile(old_key_fname):
        old_key = LoadKey(old_key_fname)
        n_old = old_key.shape[0]
    else:
        print('%s file does not exist - exiting' % old_key_fname)
        sys.exit(1)
    
    # Load new label key
    if os.path.isfile(new_key_fname):
        new_key = LoadKey(new_key_fname)
        n_new = new_key.shape[0]
    else:
        print('%s does not exist - exiting' % new_key_fname)
        sys.exit(1)

    # Check for duplicate indices and names in old and new keys        
    if CheckDuplicates(old_key, new_key):
        print('*** Exiting')
        sys.exit(1)

    # Construct label mappings
    i_old = -np.ones([n_old,])
    i_new = -np.ones([n_new,])
    
    # Init key mapping
    count = 0
    missing_key = False

    for i, old_name in enumerate(old_key.Name):
        
        old_idx = old_key.Index[i]
        new_idx = new_key.Index[new_key.Name == old_name]

        # Check where old name found in new key
        if new_idx.empty:
            print('*** %s not found in new key - skipping' % old_name)
            missing_key = True
        else:
            print('%20s: %6d -> %6d' % (old_name, old_idx, new_idx.iloc[0]))
            i_old[count] = old_idx
            i_new[count] = new_idx
            count += 1

    if missing_key:
        print('*** Could not construct a complete one-to-one label mapping - exiting')
        sys.exit(1)
            
    print('\nFound %d mappings between old and new keys' % count)
    
    # Loop over all label volumes provided
    for old_fname in label_fnames:
        
        # Construct output filename    
        old_stub, old_ext = os.path.splitext(old_fname)
        if old_ext == '.gz':
            old_stub, _ = os.path.splitext(old_stub)
        new_fname = old_stub + '_remapped.nii.gz'         
        
        print('Remapping %s to %s' % (old_fname, new_fname))
            
        # Load old label image
        old_nii = nib.load(old_fname)
        old_labels = old_nii.get_data()
        T = old_nii.get_affine()

        # Create zeroed new label image        
        new_labels = np.zeros_like(old_labels)
        
        # Loop over all old label indices
        for i, old_idx in enumerate(i_old):
            new_idx = i_new[i]
            new_labels[old_labels == old_idx] = new_idx
    
        new_nii = nib.Nifti1Image(new_labels, T)
        new_nii.to_filename(new_fname)
   
    print('Done')
    

def LoadKey(key_fname):
    '''
    Parse an ITK-SNAP label key file
    '''
    
    # Import key as a data table
    # Note the partially undocumented delim_whitespace flag
    data = pd.read_table(key_fname,
                         comment='#',
                         header=None,
                         names=['Index','R','G','B','A','Vis','Mesh','Name'],
                         delim_whitespace=True)
    
    return data


def CheckDuplicates(old_key, new_key):
    
    dups = False
    
    # Basic checks for duplicate indices and names in either key
    if old_key.duplicated('Name').any():
        dups = True
        print('*** Detected duplicate names in old key')

    if new_key.duplicated('Name').any():
        dups = True
        print('*** Detected duplicate names in new key')
    
    if old_key.duplicated('Index').any():
        dups = True
        print('*** Detected duplicate indices in new key')

    if new_key.duplicated('Index').any():
        dups = True
        print('*** Detected duplicate indices in new key')
    
    return dups

# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
