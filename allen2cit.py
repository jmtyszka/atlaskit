#!/usr/bin/python
"""
Map the Allen developmental human brain ontology (16) to CIT
simple label numbering. Output master conversion table as a CSV file.

Usage
----
allen2cit.py -o <conversion table CSV file>
allen2cit.py -h

Example
----
>>> allen2cit.py -o Allen2CIT.csv

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
import numpy as np
from urllib.request import urlopen
import xml.etree.cElementTree as ET

def main():
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Allen to CIT label index conversion')
    parser.add_argument('-o','--output', help="Conversion table CSV file")
    
    # Parse command line arguments
    args = parser.parse_args()
    
    if args.output:
        out_fname = args.output
    else:
        out_fname = 'Allen2CIT.csv'
    
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
    out_file = open(out_fname, "w")
    print('Writing conversion table to %s' % out_fname)
    
    # Write column headers
    out_file.write('CIT_Label,Allen_Label,Acronym,Name\n')
    
    for struct in tree.iter(tag='structure'):
        
            name = struct.find('name').text
            acro = struct.find('acronym').text

            # Get Allen ontology id and map to CIT simple index
            allen_id = int(struct.find('id').text)
            cit_id = np.where(ids == allen_id)[0][0]
            
            # Write line to output CSV file
            out_file.write('%d,%d,%s,%s\n' % (cit_id, allen_id, acro, name))
  

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
