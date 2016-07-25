#!/usr/bin/env python3
"""
Create an ITK-SNAP label key from a subset of the Allen Brain Atlas human
ontology via the online API

Usage
----
allen2itksnap.py -k <ITK-SNAP label key file>
allen2itksnap.py -h

Example
----
>>> allen2itksnap.py -o Amy_BG_Labels.txt

Authors
----
Mike Tyszka, Caltech, Division of Humaninities and Social Sciences

Dates
----
2015-11-20 JMT From scratch

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

import numpy as np
import argparse
from colorsys import hsv_to_rgb
from urllib.request import urlopen
import xml.etree.cElementTree as ET

def main():
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Allen ontology to ITKSNAP labels')
    parser.add_argument('-o','--output', help="Output ITKSNAP label key file")
    
    # Parse command line arguments
    args = parser.parse_args()
    
    if args.output:
        itksnap_fname = args.output
    else:
        itksnap_fname = 'Allen_Labels.txt'
    
    # Allen Brain Institute Human, 34 years, Cortex Structure Graph
    # Use developing rather than reference human brain atlas
    # http://help.brain-map.org/display/api/Atlas+Drawings+and+Ontologies
    allen_url = 'http://api.brain-map.org/api/v2/structure_graph_download/16.xml'
    
    # Download XML for Allen human brain ontology
    print('Downloading ontology from %s' % allen_url)
    xml_remote = urlopen(allen_url)

    # Parse XML as a Document Object Model (DOM) 
    print('Parsing XML tree')
    tree = ET.parse(xml_remote)
    
    # Collect and sort all label numbers into an array
    ids = list()
    for id in tree.iter(tag='id'):
        ids.append(int(id.text))
    ids = np.sort(np.array(ids))
    
    # Open output text file
    out_file = open(itksnap_fname, "w")
    print('Writing label key to %s' % itksnap_fname)
    
    # List of structures to expand
    # - Cerebral Nuclei (CN). Amygdala, etc
    # - Diencephalon (Die). Thalamus, etc
    # - Forebrain White Matter (FWM). AC, etc
    # - Midbrain Tegmentum (MTg). SN, etc
    acronyms = ['"CN"','"Die"','"FWM"','"MTg"']
    
    for ac, acro in enumerate(acronyms):
        
        # Cycle through hues for each structure group
        # Color-anticolor pairs via Ha
        # Generates H = 0, 180, 60, 240, 120, 300 then repeats
        Ha = np.mod(ac,2)
        H = np.mod(int(ac/2) * 60.0 + Ha * 180.0, 360.0) / 360.0
        
        print('\nExpanding %s :' % acro)
        
        for struct in tree.iter(tag='structure'):
            
            if struct.find('acronym').text == acro:
                
                # Found structure with correct acronym
                # Iterate into this structure, printing out info
                
                for sc, substruct in enumerate(struct.iter(tag='structure')):
                    
                    # Structure acronym from Allen ontology
                    sub_acro = substruct.find('acronym').text
                    
                    print(sub_acro, end=' ')

                    # Get Allen id and map to CIT id
                    allen_id = int(substruct.find('id').text)
                    cit_id = np.where(ids == allen_id)[0][0]
                    
                    # Convert CIT id to RGB triple via HSV

                    # S = 0.25:0.10:0.75 [6 levels] (outer loop)
                    S = np.mod(cit_id/6, 6) * 0.1 + 0.25

                    # V = 0.50:0.10:1.00 [6 levels] (inner loop)
                    V = np.mod(cit_id, 6) * 0.1 + 0.50
                    
                    # Convert HSV to RGB triple
                    R, G, B = hsv_to_rgb(H, S, V)
                    
                    # Scale RGB to [0,255] range
                    R, G, B = int(R*255.0), int(G*255.0), int(B*255.0)
                                        
                    # Write line to output ITKSNAP label file
                    # print('%5d%6d%6d%6d%9d%3d%3d     %s' % (cit_id, R, G, B, 1, 1, 1, sub_acro))
                    out_file.write('%5d%6d%6d%6d%9d%3d%3d     %s\n' % (cit_id, R, G, B, 1, 1, 1, sub_acro))
                
                # End substructure acronym list with newline
                print('')

    # Clean up
    out_file.close()
    print('Done')


def Hex2RGB(sHex):
    '''
    Convert RGB hexadecimal triplet string into separate decimal values
    '''
    
    # Divide into substrings and convert to decimal integers
    R = int(sHex[0:2], 16)
    G = int(sHex[2:4], 16)
    B = int(sHex[4:6], 16)
    
    return R, G, B


def SaveKey(key_fname, ontology):
    '''
    Save an ontology as an ITK-SNAP label key file
    '''
    
    print('Saving ontology in %s' % key_fname)



# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
