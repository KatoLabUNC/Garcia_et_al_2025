% Convert structure.nrrd file to surface mesh data.
%  Hiroyuki Kato. 
%
% Matlab 2022b and later is needed.
% Before running the code, you need to download structure_nrrd files corresponding to IC subdivisions from Allen CCF.
% Download structure mask files structure_4.nrrd, structure_811.nrrd, structure_820.nrrd, and structure_828.nrrd, which correspond to IC, ICc, ICd, and ICe, respectively,
% from https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/ccf_2017/structure_masks/structure_masks_10/ .
%

% Define the path to your Allen CCF directory
ccfPath = '(path where you saved the nrrd files)';

% pick which mesh to make
structure_to_map = 'IC';
structure_to_map = 'ICc';
structure_to_map = 'ICd';
structure_to_map = 'ICe';


% Full path to the .nrrd file
switch structure_to_map
    case 'IC'
        nrrdFilePath = fullfile(ccfPath, 'structure_masks_10\', 'structure_4.nrrd');
    case 'ICc'
        nrrdFilePath = fullfile(ccfPath, 'structure_masks_10\', 'structure_811.nrrd');
    case 'ICd'
        nrrdFilePath = fullfile(ccfPath, 'structure_masks_10\', 'structure_820.nrrd');
    case 'ICe'
        nrrdFilePath = fullfile(ccfPath, 'structure_masks_10\', 'structure_828.nrrd');
end

% Read the .nrrd file using the built-in nrrdread function
data = nrrdread(nrrdFilePath); % This is a mask corresponding to structure_4, which is IC. 1320 x 800 x 1140 (AP x DV x ML). Note, in Allen ccf, AP is x, DV is y, ML is z.

% Verify that data has been loaded
if isempty(data)
    error('Failed to read data from the .nrrd file.');
end

% Check the unique values in the data
uniqueValues = unique(data(:));
disp('Unique data values in the NRRD file:');
disp(uniqueValues');

% If the data is binary (e.g., contains 0 and 1), we can proceed to create a mask
% If not, adjust the threshold accordingly

% Assuming the data contains binary values where the structure is represented by 1
% Create a binary mask for the inferior colliculus
mask = data > 0;

% Check if the mask contains any true values
if ~any(mask(:))
    error('No voxels found in the data.');
end

% show only right hemisphere
mask(:,:,1:570) = 0;


% Visualize the mask using isosurface
[faces, vertices] = isosurface(mask, 0.5);
faces = faces(:,[2 3 1]);
vertices = vertices(:,[2 3 1]); % somehow, isosurface switches axes

% Plot the isosurface
figure;
patch('Vertices', vertices, 'Faces', faces, ...
      'FaceColor', 'cyan', 'EdgeColor', 'none', 'FaceAlpha', 0.3);

% Enhance visualization
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Inferior Colliculus Mesh from NRRD File');
view(3);
camlight('headlight');
lighting gouraud;
material dull;
grid on;


% %% create a mask for IC
% 
% nrrdFilePath = fullfile(ccfPath, 'annotation_10.nrrd');
% 
% % Read the .nrrd file using the built-in nrrdread function
% data = nrrdread(nrrdFilePath);
% mask = data == 4; % 4 is IC's structure ID
% 
% save('(save path)\IC_mask.mat','mask');

%% save 3D mesh


switch structure_to_map
    case 'IC'
        save('(save path)\Allen_ccf\structure_masks_10\IC_mesh.mat','faces','vertices');
    case 'ICc'
        save('(save path)\Allen_ccf\structure_masks_10\ICc_mesh.mat','faces','vertices');
    case 'ICd'
        save('(save path)\Allen_ccf\structure_masks_10\ICd_mesh.mat','faces','vertices');
    case 'ICe'
        save('(save path)\Allen_ccf\structure_masks_10\ICe_mesh.mat','faces','vertices');
end
