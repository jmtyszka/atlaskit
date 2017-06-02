#!/usr/bin/env python3
"""
Download the entire Allen Reference Human Brain Atlas (34 years, cortex-gyral) as a JPG stack
- use the brain-map.org API

Usage
----
allen_download.py

Authors
----
Mike Tyszka, Caltech, Division of Humaninities and Social Sciences

Dates
----
2017-05-24 JMT From scratch

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
2017 California Institute of Technology.
"""

__version__ = '0.1.0'

import os
import sys
import requests
import xml.etree.ElementTree as ET
from urllib.request import urlopen


def main():

    # Image type: 'vector' (SVG) or 'image' (JPEG, annotated)
    imtype = 'image'
    
    # Allen Brain Institute Human, 34 years, Cortex-Gyral atlas
    # Use developing rather than reference human brain atlas
    # http://help.brain-map.org/display/api/Atlas+Drawings+and+Ontologies

    # API request for full list of all images associated with the Human 34y gyral atlas
    api_request_url = ("http://api.brain-map.org/api/v2/data/query.xml?criteria=model::AtlasImage,"
    "rma::criteria,[annotated$eqtrue],"
    "atlas_data_set(atlases[id$eq138322605]),"
    "alternate_images[image_type$eq'Atlas+-+Developing+Human'],"
    "rma::options[order$eq'sub_images.section_number'][num_rows$eqall]")

    # Base image download format string
    if imtype == 'vector':
        image_url_fmt = 'http://api.brain-map.org/api/v2/svg_download/%s?groups=31,113753815&downsample=8d'
        out_dir = '/Users/jmt/sandbox/Allen_Vector'
    else:
        image_url_fmt = 'http://api.brain-map.org/api/v2/atlas_image_download/%s?groups=31,113753815&downsample=4&annotation=true'
        out_dir = '/Users/jmt/sandbox/Allen_Image'

    # Safely create output directory
    os.makedirs(out_dir, exist_ok=True)
    
    # Download XML for Allen human brain ontology
    print('Downloading XML image list from brain-map.org')
    try:
        image_list_xml = urlopen(api_request_url)
    except:
        print('* Problem getting XML image list from server')
        sys.exit(1)

    # Parse XML as a Document Object Model (DOM) 
    print('Parsing XML image list')
    tree = ET.parse(image_list_xml)
    
    # Loop over all atlas images
    for atlas_image in tree.iter(tag='atlas-image'):

        id_tag = atlas_image.find('id')
        sno_tag = atlas_image.find('section-number')
        sno = int(sno_tag.text)

        # Construct image download URL
        image_url = image_url_fmt % id_tag.text

        # Construct output filename
        if imtype == 'vector':
            image_fname = os.path.join(out_dir, 'Allen_%04d.svg' % sno)
        else:
            image_fname = os.path.join(out_dir, 'Allen_%04d.jpg' % sno)

        print('  Writing image to %s' % image_fname)

        # Grab image data and write to a files
        image_data = requests.get(image_url).content
        with open(image_fname, 'wb') as handler:
            handler.write(image_data)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
