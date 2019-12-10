function fs_roi_overlay(roi_csv)
% FS_ROI_OVERLAY Create a surface overlay from ROI results
%   Map scalar values for each label ROI in the Destrieux atlas
%   to a Freesurfer curvature overlay for use by Freeview
%
% AUTHOR : Mike Tyszka
% PLACE  : Caltech
% DATES  : 2019-07-16 JMT From scratch
%          2019-07-18 JMT Switch to curv file output

if nargin < 1
  fprintf('USAGE : fs_roi_surf(roi_csv)\n');
  return
end

% Construct output annotations filename
[fpath, fstub] = fileparts(roi_csv);
lh_curv_fname = fullfile(fpath, [fstub '_lh.curv']);
rh_curv_fname = fullfile(fpath, [fstub '_rh.curv']);

% Find Freesurfer Average directory
fsavg_dir = fullfile('/Applications/freesurfer/subjects/fsaverage');
fsavg_label_dir = fullfile(fsavg_dir, 'label');

% Destrieux 2009 hemisphere annotation files
fsavg_lh_annot = fullfile(fsavg_label_dir, 'lh.aparc.a2009s.annot');
fsavg_rh_annot = fullfile(fsavg_label_dir, 'rh.aparc.a2009s.annot');

% Load the Destrieux 2009 annotation data
% Vertices : vertex index into fsaverage surface
% Labels : integer parcel labels
% Color Table : table mapping integer parcel labels to parcel name & RGBA
[lh_verts, lh_labels, lh_ctable] = read_annotation(fsavg_lh_annot);
[rh_verts, rh_labels, rh_ctable] = read_annotation(fsavg_rh_annot);

% Load the label names and values from the csv file
[lh_roi_names, lh_roi_vals, rh_roi_names, rh_roi_vals] = load_roi_data(roi_csv);

% Map ROI values to vertices in fsaverage surface
lh_curv = map_to_surface(lh_roi_names, lh_roi_vals, lh_verts, lh_labels, lh_ctable);
rh_curv = map_to_surface(rh_roi_names, rh_roi_vals, rh_verts, rh_labels, rh_ctable);

% Write curvature-format files for LH and RH
write_curv(lh_curv_fname, lh_curv, length(lh_curv));
write_curv(rh_curv_fname, rh_curv, length(rh_curv));


function vals = map_to_surface(roi_names, roi_vals, verts, labels, ctable)
% Map ROI scalar values to atlas label regions

vals = zeros(size(verts));

% Loop over ROI names
for rc = 1:length(roi_names)

  % Find ROI name in color table
  mask = ismember(ctable.struct_names, roi_names{rc});
  
  if sum(mask) > 0
  
    % Pull label index from colortable
    label_idx = ctable.table(mask, 5);

    % Find all vertices in annotation with this label index
    vert_mask = (labels == label_idx);
  
    % Set all vertices within label to ROI scalar value
    vals(vert_mask) = roi_vals(rc);
    
  end
  
end

return


function [lh_names, lh_vals, rh_names, rh_vals] = load_roi_data(roi_csv)

T = readtable(roi_csv, 'delimiter', ',');

roi_names = T{:, 1};
roi_vals = T{:, 2};

% Split into LH and RH name-val pairs
lh_mask = contains(roi_names, 'Left');
rh_mask = contains(roi_names, 'Right');

lh_names = preproc_names(roi_names(lh_mask));
lh_vals = roi_vals(lh_mask);

rh_names = preproc_names(roi_names(rh_mask));
rh_vals = roi_vals(rh_mask);

return



function new_names = preproc_names(old_names)
% Convert to ctable stuct name format for searches

% Remove the hemisphere prefix
new_names = erase(old_names, 'Right.');
new_names = erase(new_names, 'Left.');

% Replace 'G.S' with 'G_and_S'
new_names = replace(new_names, 'G.S', 'G_and_S');

% Replace remaining '.'s with '-'s
new_names = replace(new_names, '.', '-');

return