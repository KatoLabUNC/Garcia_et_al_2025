%% Create dataset of brain section images.
% Koun Onodera and Hiroyuki Kato. 220615
%
% This code creates SectionImage.mat.
% Before running this code, you need to prepare tif files including the region of interest within the brain atlas.
% For our MGB analysis, we used Paxinos Mouse Brain Atlas, and included the following range:
% AP: Bregma -4.04mm through -2.80mm
% ML: Lateral 1.50mm through 2.50mm
% DV: Ventral 2.50mm through 4.00mm.
%
% For example, Bregma-2.80mm_BrainSection.tif shows a traced brain atlas at Bregma -2.80mm, in the range of lateral 1.5-2.5mm, ventral 2.5-4mm.
% If you have access to the digital version of the atlas, simply use cliping masks to extract that portion of the atlas data. 
%
% Note that we cannot share these files on GitHub due to copyright restrictions.
% Please prepare these files yourself or contact Hiroyuki Kato.
%
% SectionImage is a 1x11 struct with 4 fields (APposition, Image, Xlength, and Ylength).
% 11 entries correspond to 11 atlas brain sections containing MGB (Bregma -2.80 through -4.04).
% APposition: in mm. (e.g. -2.8000)
% Image: brain atlas image including MGB. (e.g., 893x594 uint8 matrix, made by tracing the Paxinos Atlas. In our case, we traced the area within X = 1.5-2.5mm, Y = 2.5-4mm.)
% Xrange: range of X included in Image. (e.g., [1.5000, 2.5000] would be sufficient to contain the entire MGB.)
% Yrange: range of Y included in Image. (e.g., [2.5000, 4.0000] would be sufficient to contain the entire MGB.);

clear

filepath = '(filepath that stores tif files of traced brain atlas section images)';


%% load MGB images and save as a mat file

APpos = [-2.80 -2.92 -3.08 -3.16 -3.28 -3.40 -3.52 -3.64 -3.80 -3.88 -4.04]; % These are AP positions of Paxinos Brain Atlas sections containing MGB.

for i = 1:length(APpos)
    SectionImage(i).APposition = APpos(i);
    SectionImage(i).Image = imread([filepath filesep 'Bregma' num2str(APpos(i),'%.2f') 'mm_BrainSection.tif']);
    SectionImage(i).Xrange = [1.5 2.5];
    SectionImage(i).Yrange = [2.5 4];
end

save([filepath filesep 'SectionImage.mat'], 'SectionImage');
