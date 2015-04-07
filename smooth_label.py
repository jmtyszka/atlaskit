#!/opt/local/Library/Frameworks/Python.framework/Versions/2.7/Resources/Python.app/Contents/MacOS/Python
"""
Smooth a single label within an integer atlas volume

Usage
----
smooth_label.py <source atlas> <smoothed atlas> <label number>

Example
----
>>> smooth_label.py atlas.nii.gz atlas_smooth_5.nii.gz 5

Authors
----
Mike Tyszka, Caltech Brain Imaging Center

Dates
----
2015-04-07 JMT From scratch

License
----
This file is part of mrgaze.

    atlaskit is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    atlaskit is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with mrgaze.  If not, see <http://www.gnu.org/licenses/>.

Copyright
----
2015 California Institute of Technology.
"""

__version__ = '0.0.1'

import os
import sys


def main():
    
    # Get command line arguments
    if len(sys.argv) < 3:
        print('Usage: smooth_label.py <source atlas> <smoothed atlas> <label number>')
        sys.exit(1)
    else:
        src_img = sys.argv[1]
        smth_img = sys.argv[2]
        label = sys.argv[3]
        
    # Load the source atlas image
    src_nii = 
    
    # Clean exit
    sys.exit(0)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()