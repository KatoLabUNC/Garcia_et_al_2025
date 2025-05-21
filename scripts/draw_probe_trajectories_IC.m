% draw probe trajectories of linear probe recordings.
% Hiroyuki Kato, 240218
%
% Before running this code, run AP_histology first. 
%

clear all


%% define paths
parentdir = '(parent directory path)';



% inferior colliculus recordings 
m = 1; 
mouse = []; date = []; site = []; target_ch = []; condition = [];

mouse{m} = '(mouse1)'; date{m} = '(date1)'; site{m} = '(site1)'; m = m+1;
mouse{m} = '(mouse2)'; date{m} = '(date2)'; site{m} = '(site2)'; m = m+1;
mouse{m} = '(mouse3)'; date{m} = '(date3)'; site{m} = '(site3)'; m = m+1;
% add as many experiments as you want


%% load atlas

allen_atlas_path = fileparts(which('template_volume_10um.npy')); % This file comes with AP_histology package.
if isempty(allen_atlas_path)
    error('No CCF atlas found (add CCF atlas to path)')
end
disp('Loading Allen CCF atlas...')
tv = readNPY(fullfile(allen_atlas_path,'template_volume_10um.npy'));
av = readNPY(fullfile(allen_atlas_path,'annotation_volume_10um_by_index.npy'));
st = ap_histology.loadStructureTree(fullfile(allen_atlas_path,'structure_tree_safe_2017.csv'));
disp('Done.')

%% go through mice

figure('Name','Probe trajectories');
axes_atlas = axes;
[~, brain_outline] = plotBrainGrid([],axes_atlas); % change line color in plotBrainGrid as necessary. plotBrainGrid comes with AP_histology package.
[ap_max,dv_max,ml_max] = size(tv);
xlim([-10,ap_max+10])
ylim([-10,ml_max+10])
zlim([-10,dv_max+10])

h = rotate3d(gca);
h.Enable = 'on';


% plot IC
load('(your file path for IC mesh files)\Allen_ccf\structure_masks_10\ICc_mesh.mat'); 
patch('Vertices', vertices, 'Faces', faces, ...
      'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.4);
load('(your file path for IC mesh files)\Allen_ccf\structure_masks_10\ICd_mesh.mat'); 
patch('Vertices', vertices, 'Faces', faces, ...
      'FaceColor', 'green', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
load('(your file path for IC mesh files)\Allen_ccf\structure_masks_10\ICe_mesh.mat');
patch('Vertices', vertices, 'Faces', faces, ...
      'FaceColor', 'cyan', 'EdgeColor', 'none', 'FaceAlpha', 0.2);


for m = 1:length(mouse)
    
    subdir = [parentdir filesep date{m} filesep mouse{m} filesep 'processed_img']; % our folder structure is parentdir/date/mouse/processed_img/probe_ccf.mat
        
    % Load corresponding probe ccf file
    load(fullfile(subdir,'probe_ccf.mat'));
    
    
    % Plot line of best fit
    r0 = mean(probe_ccf.points,1);
    xyz = bsxfun(@minus,probe_ccf.points,r0);
    [~,~,V] = svd(xyz,0);
    histology_probe_direction = V(:,1);
    % (make sure the direction goes down in DV - flip if it's going up)
    if histology_probe_direction(2) < 0
        histology_probe_direction = -histology_probe_direction;
    end

    line_eval = [-1000,1000];
    probe_fit_line = bsxfun(@plus,bsxfun(@times,line_eval',histology_probe_direction'),r0);
    
    
    % if left hemisphere, flip to the right hemisphere
    % Here, 1st dim: AP, 2nd dim: DV, 3rd dim: ML. ML = 570 corresponds to midline.
    if mean(probe_fit_line(:,3))<570
        probe_fit_line(:,3) = 570-(probe_fit_line(:,3)-570);
    end
    
    
    % optionally, ICc in red , ICe in blue, and ICd in green
    
    % restrict the line only between marked points. HK 240218
    % use Z value to find edges of distribution
    Z_min = min(probe_ccf.points(:,2));
    Z_max = max(probe_ccf.points(:,2));
    depth_min = probe_ccf.trajectory_areas.trajectory_depth(1,1);
    depth_max = probe_ccf.trajectory_areas.trajectory_depth(end,2);

    non_IC_segments = find(cell2mat(cellfun(@(x) isempty(strfind(x,'ICc'))...
                                                       && isempty(strfind(x,'ICe'))...
                                                       && isempty(strfind(x,'ICd')), probe_ccf.trajectory_areas.acronym,'uniformoutput',0)));
   

    ICc_segments = find(cell2mat(cellfun(@(x) ~isempty(strfind(x,'ICc')),probe_ccf.trajectory_areas.acronym,'uniformoutput',0)));
    
    for s = 1:length(ICc_segments)
        curr_depth = probe_ccf.trajectory_areas.trajectory_depth(ICc_segments(s),:);
        curr_Z = interp1([depth_min depth_max],[Z_min Z_max],curr_depth);
        probe_fit_segment_interp = [interp1(probe_fit_line(:,2),probe_fit_line(:,1),curr_Z)' ...
                                    curr_Z' ...
                                    interp1(probe_fit_line(:,2),probe_fit_line(:,3),curr_Z)'];
        line(probe_fit_segment_interp(:,1),probe_fit_segment_interp(:,3),probe_fit_segment_interp(:,2), ...
             'color','r','linewidth',1.5)
    end
    
    ICe_segments = find(cell2mat(cellfun(@(x) ~isempty(strfind(x,'ICe')),probe_ccf.trajectory_areas.acronym,'uniformoutput',0)));
    
    for s = 1:length(ICe_segments)
        curr_depth = probe_ccf.trajectory_areas.trajectory_depth(ICe_segments(s),:);
        curr_Z = interp1([depth_min depth_max],[Z_min Z_max],curr_depth);
        probe_fit_segment_interp = [interp1(probe_fit_line(:,2),probe_fit_line(:,1),curr_Z)' ...
                                    curr_Z' ...
                                    interp1(probe_fit_line(:,2),probe_fit_line(:,3),curr_Z)'];
        line(probe_fit_segment_interp(:,1),probe_fit_segment_interp(:,3),probe_fit_segment_interp(:,2), ...
             'color',[0 0.3 0.7],'linewidth',1.5)
    end
    
    ICd_segments = find(cell2mat(cellfun(@(x) ~isempty(strfind(x,'ICd')),probe_ccf.trajectory_areas.acronym,'uniformoutput',0)));
    
    for s = 1:length(ICd_segments)
        curr_depth = probe_ccf.trajectory_areas.trajectory_depth(ICd_segments(s),:);
        curr_Z = interp1([depth_min depth_max],[Z_min Z_max],curr_depth);
        probe_fit_segment_interp = [interp1(probe_fit_line(:,2),probe_fit_line(:,1),curr_Z)' ...
                                    curr_Z' ...
                                    interp1(probe_fit_line(:,2),probe_fit_line(:,3),curr_Z)'];
        line(probe_fit_segment_interp(:,1),probe_fit_segment_interp(:,3),probe_fit_segment_interp(:,2), ...
             'color',[0 0.7 0.3],'linewidth',1.5)
    end
    
    % note, 1st dim: AP, 2nd dim: ML, 3rd dim: DV. ML = 570 corresponds to midline.
    
end

xlim([900 1230]);
ylim([570 900]);
zlim([50 380]);

% for saving data, use
view([-150 30])
% view([180 0])
view([270 0])
% view([90 0])
view([0 90])
% and save by 
% export_fig (filename) -transparent
