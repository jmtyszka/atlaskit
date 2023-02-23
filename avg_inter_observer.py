#!/usr/bin/env python3
"""
Create report of average intra and inter-observer atlas label statistics

Usage
----
avg_inter_observer.py 


Authors
----
Wolfgang M. Pauli, Caltech Brain Imaging Center

Dates
----
2017-10-27 WMP Split from atlas.py

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

from atlas_report import *

atlas_dir = os.environ['ATLAS_DIR'] # directory containing e.g. inter_observer_metrics.csv

intra_stats, inter_stats = load_metrics(atlas_dir)
report_dir = '/tmp/'
label_names, label_nos, observers, templates, dice, haus = inter_stats

ncols = 4
nrows = np.ceil(len(label_names)/ncols).astype(int)

# Metric limits
dlims = 0.0, 1.0
hlims = 0.0, 10.0
hlims_sd = 0.0, 2.0
dlims_sd = 0.0, .2

# Init image filename lists for HTML template
inter_dice_imgs = []
inter_haus_imgs = []

# Create similarity figures over all labels and observers
dice_fname = "inter_tmp_dice_m.png" 
dice_m = dice.mean(1)
similarity_figure(dice_m[:,:,:],
                  "Dice Coefficient (Means)",
                  dice_fname,
                  report_dir, label_names, dlims, nrows, ncols, 0.0, 12, True)

    
haus_fname = "inter_tmp_haus_m.png"
haus_m = haus.mean(1)
similarity_figure(haus_m[:,:,:],
                  "Hausdorff Distance (Means)",
                  haus_fname,
                  report_dir, label_names, hlims, nrows, ncols, 1e6, 12, True)

    
# Create similarity figures over all labels and observers
dice_fname = "inter_tmp_dice_sd.png" 
dice_sd = np.sqrt(dice.var(1))
similarity_figure(dice_sd[:,:,:],
                  "Dice Coefficient (SD)",
                  dice_fname,
                  report_dir, label_names, dlims_sd, nrows, ncols, 0.0, 12, True)

    
haus_fname = "inter_tmp_haus_sd.png"
haus_sd = np.sqrt(haus.var(1))
similarity_figure(haus_sd[:,:,:],
                  "Hausdorff Distance (SD)",
                  haus_fname,
                  report_dir, label_names, hlims_sd, nrows, ncols, 1e6, 12, True)
