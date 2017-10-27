from atlas_report import *

atlas_dir = os.environ['ATLAS_DIR']

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
                  report_dir, label_names, dlims, nrows, ncols, 0.0, 10)

    
haus_fname = "inter_tmp_haus_m.png"
haus_m = haus.mean(1)
similarity_figure(haus_m[:,:,:],
                  "Hausdorff Distance (Means)",
                  haus_fname,
                  report_dir, label_names, hlims, nrows, ncols, 1e6, 10)

    
# Create similarity figures over all labels and observers
dice_fname = "inter_tmp_dice_sd.png" 
dice_sd = np.sqrt(dice.var(1))
similarity_figure(dice_sd[:,:,:],
                  "Dice Coefficient (SD)",
                  dice_fname,
                  report_dir, label_names, dlims_sd, nrows, ncols, 0.0, 10)

    
haus_fname = "inter_tmp_haus_sd.png"
haus_sd = np.sqrt(haus.var(1))
similarity_figure(haus_sd[:,:,:],
                  "Hausdorff Distance (SD)",
                  haus_fname,
                  report_dir, label_names, hlims_sd, nrows, ncols, 1e6, 10)
