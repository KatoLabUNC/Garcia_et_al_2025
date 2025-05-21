% save layer info using probe_ccf file
% Hiroyuki Kato. 240703
%
% Before running this code, run AP_histology (https://github.com/petersaj/AP_histology) to align probe trajectories to the Allen Common Coordinate Framework.
% This code generates probe_ccf.mat, which we use to assign each linear probe channel to brain regions.
%
% Additionally, before running this code, identify the easily-identifiable landmark channels within your linear probe recording.
% This would include the brain surface, the boundary between the cortex and the inferior colliculus, etc. 
% Next, within the probe_ccf.mat, go to probe_ccf.trajectory_areas, and find the depths corresponding to these landmarks under "trajectory_depth".
% Provide the trajectory depths as the first column, and the landmark channels as the second column of "anchorpoints".
%


%% define paths and parameters

parentdir = '(parent directory)';

mouse = '(mouse name)'; % data folder is organized as parentdir\date\mouse\site\
date = '(experiment date)';
site = '(recording site ID)';
surfaceCh = 5; % example. channel corresponding to the surface
totalCh = 286; % example. number of recorded channels
anchorpoints = [1244 124; 1412 150]; % example. see the explanation above.


% anchorpoints show ccf depths on the 1st column and channel number on the 2nd column.
% Anchor points can be things like the cortex-midbrain border, cortex-white matter border, etc.,
% which can be determined by spike distributions.


%% atlas area names
                 
% load ccf file
load([parentdir filesep date filesep mouse filesep 'processed_img' filesep 'probe_ccf.mat']); % move probe_ccf.mat to appropriate directory.

trajectory_origin = probe_ccf.trajectory_coords(1,:)*10; % This is the CCF coordinates of the probe entry point. converted to um.
% note that coordinates are (AP, DV, ML), and ML=5700 corresponds to the midline. (0-5700: left hemisphere. 5700-11400: right hemisphere).
% The origin is rostral, dorsal, left corner of the 3d volume, and bregma is roughly (5400, 440, 5700).

trajectory_vector = probe_ccf.trajectory_coords(2,:)-probe_ccf.trajectory_coords(1,:); 
trajectory_vector = trajectory_vector/norm(trajectory_vector); % vector with length 1.Use this, together with the depth information, to estimate the coordinates.

probe_ccf = probe_ccf.trajectory_areas;

trajectory_borders = probe_ccf.trajectory_depth;
trajectory_borders = [trajectory_borders(:,1); trajectory_borders(end,2)];

trajectory_areas = probe_ccf.acronym;


%% find right subfolder


subdir = [parentdir filesep date filesep mouse filesep site];

layer_info = [];
layer_info.mouse = mouse;
layer_info.date = date;
layer_info.site = site;
layer_info.surfaceCh = surfaceCh;
layer_info.totalCh = totalCh;

anchorpoints = [0 surfaceCh;
                anchorpoints;
                trajectory_borders(end) totalCh];

            
% calculate coordinates of each anchor points.
anchor_coordinates = cell(length(anchorpoints),1);
for anch = 1:size(anchor_coordinates,1)
    anchor_coordinates{anch} = trajectory_origin + anchorpoints(anch,1)*trajectory_vector;
end

layer_info.ch_label = cell(totalCh,1);
% go through channels between fixed anchorpoints
for anch = 1:size(anchorpoints,1)-1
    
    curr_edges = [anchorpoints(anch,1) anchorpoints(anch+1,1)];
    curr_channels = anchorpoints(anch,2):anchorpoints(anch+1,2);
    if isempty(curr_channels)
        continue
    end
    
    % if the anchorpoint is the brain surface, add 0.5 to the channel number, as our definition of surface
    % is just outside the brain.
    if anch == 1
        channel_depths = interp1([curr_channels(1)+0.5 curr_channels(end)], curr_edges, curr_channels(2:end));
        curr_channels(1) = [];
    else
        channel_depths = interp1([curr_channels(1) curr_channels(end)], curr_edges, curr_channels);
    end
    
    for ch = 1:length(curr_channels)
        curr_depth = channel_depths(ch);
        segment_id = find(curr_depth> trajectory_borders, 1, 'last');
        % make sure that its within the segment
        segment_id2 = find(curr_depth<=trajectory_borders, 1, 'first')-1;
        if segment_id~=segment_id2
            warning(['corresponding segment not found for ch' num2str(curr_channels(ch))]);
            continue
        end
        
        layer_info.ch_label{curr_channels(ch)} = trajectory_areas{segment_id};
        
        % coordinates in CCF
        layer_info.ch_coordinates_CCF{curr_channels(ch)} = anchor_coordinates{1}+curr_depth*trajectory_vector;
        
        % convert to bregma-based coordinates
        coord_from_bregma = layer_info.ch_coordinates_CCF{curr_channels(ch)}-[5400 440 5700];
        coord_from_bregma(1) = -coord_from_bregma(1);
        layer_info.ch_coordinates_bregma{curr_channels(ch)} = coord_from_bregma;
    end
end

save([subdir filesep 'layer_info2.mat'],'layer_info');

% layer_info.mat contains the brain region identity and the stereotaxic coordinate for each channel. 
% This file will be loaded in all subsequent analyses to determine properties of different IC subregions.
