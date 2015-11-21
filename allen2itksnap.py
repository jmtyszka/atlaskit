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

import sys
import argparse
from bs4 import BeautifulSoup as bs

def main():
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Create ITK-SNAP label key from Allen Human Brain ontology')
    parser.add_argument('-k','--key', help="ITK-SNAP label key file")
    
    # Parse command line arguments
    args = parser.parse_args()
    key_fname = args.key
    
    # Allen Institute human brain ontology (10)
    allen_url = 'http://api.brain-map.org/api/v2/structure_graph_download/10.xml'
    
    soup = bs(allen_url)
    
    print(soup.find_all('amygdala'))
    


def SaveKey(key_fname, ontology):
    '''
    Save an ontology as an ITK-SNAP label key file
    '''
    
    print('Saving ontology in %s' % key_fname)



# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
