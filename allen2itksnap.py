#!/usr/bin/python
"""
Create an ITK-SNAP label key from a subset of the Allen Brain Atlas human
ontology availble through the online API

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

import argparse
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
    
    # Open output text file
    out_file = open(itksnap_fname, "w")
    print('Writing label key to %s' % itksnap_fname)
    
    # List of structures to expand
    acronyms = ['"BN"','"MTg"','"AMY"','"ATA"','"FWM"']
    
    for acro in acronyms:
    
        for struct in tree.iter(tag='structure'):
            
            if struct.find('acronym').text == acro:
                
                # Found structure with correct acronym
                # Iterate into this structure, printing out info
                
                for substruct in struct.iter(tag='structure'):
                    
                    # sub_acro = substruct.find('acronym').text
                    sub_name = substruct.find('name').text
                    sub_id = int(substruct.find('id').text)
                    sub_rgb = substruct.find('color-hex-triplet').text
                    
                    # Split hex triplet into decimal R, G and B
                    R, G, B = Hex2RGB(sub_rgb)
                    
                    # Write line to output ITKSNAP label file
                    if sub_id > 99999:
                        print('%5d%6d%6d%6d%9d%3d%3d     %s' % (sub_id, R, G, B, 1, 1, 1, sub_name))
                    else:
                        out_file.write('%5d%6d%6d%6d%9d%3d%3d     %s\n' % (sub_id, R, G, B, 1, 1, 1, sub_name))
              

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
