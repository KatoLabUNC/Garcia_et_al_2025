%% Extract tissue outlines for individual MGB regions.
% Koun Onodera and Hiroyuki Kato. 220615
%
% This code generates polygons showing MGB subdivision boundaries. 
%
% Before running this code, you need to prepare tif files showing MGB subregions within the brain atlas.
% For our MGB analysis, we used Paxinos Mouse Brain Atlas, and included the following range:
% AP: Bregma -4.04mm through -2.80mm
% ML: Lateral 1.50mm through 2.50mm
% DV: Ventral 2.50mm through 4.00mm.
%
% For example, Bregma-2.80mm_BrainSection.tif shows a traced MGB boundaries of the brain atlas at Bregma -2.80mm, in the range of lateral 1.5-2.5mm, ventral 2.5-4mm.
% If you have access to the digital version of the atlas, simply use cliping masks to extract that portion of the atlas data and trace MGB subdomains. 
%
% Note that we cannot share these files on GitHub due to copyright restrictions.
% Please prepare these files yourself or contact Hiroyuki Kato.
%
% MGBdrawdata is a 1x11 struct with 7 fields (APposition, MGv, MGd, MGm, SG, BIC, and MZMG).
% 11 entries correspond to 11 atlas brain sections containing MGB.
% APposition: in mm. (e.g. -2.8000)
% MGv, MGd, MGm, SG, BIC, and MZMG fields include polygon coordinates of individual MGB subdivisions (represented in mm) within the corresponding section.
% For example, AP=-2.8000 would have polygons for MGv, MGd, MGm, and SG, while AP = -3.80000 would have polygons for BIC and MZMG.
%

clear all

filepath = '(filepath that stores tif files of traced brain atlas section images)';

%% load and show MGB images

APpos = [-2.80 -2.92 -3.08 -3.16 -3.28 -3.40 -3.52 -3.64 -3.80 -3.88 -4.04]; % These are AP positions of Paxinos Brain Atlas sections containing MGB.
xrange = [1.5 2.5];
yrange = [4 2.5];

RegionColor{1} = [1 0 0];%red
RegionColor{2} = [0 1 0];%green
RegionColor{3} = [0 0 1];%blue
RegionColor{4} = [1 0 1];%magenta
RegionColor{5} = [0 1 1];%cyan

MGBimage{1} = imread([filepath filesep 'Bregma-2.80mm_MGB.tif']);
MGBimage{2} = imread([filepath filesep 'Bregma-2.92mm_MGB.tif']);
MGBimage{3} = imread([filepath filesep 'Bregma-3.08mm_MGB.tif']);
MGBimage{4} = imread([filepath filesep 'Bregma-3.16mm_MGB.tif']);
MGBimage{5} = imread([filepath filesep 'Bregma-3.28mm_MGB.tif']);
MGBimage{6} = imread([filepath filesep 'Bregma-3.40mm_MGB.tif']);
MGBimage{7} = imread([filepath filesep 'Bregma-3.52mm_MGB.tif']);
MGBimage{8} = imread([filepath filesep 'Bregma-3.64mm_MGB.tif']);
MGBimage{9} = imread([filepath filesep 'Bregma-3.80mm_MGB.tif']);
MGBimage{10} = imread([filepath filesep 'Bregma-3.88mm_MGB.tif']);
MGBimage{11} = imread([filepath filesep 'Bregma-4.04mm_MGB.tif']);

figure('position',[50 200 1600 600]);
for i = 1:size(MGBimage,2)
    subplot(2,6,i)
    image(flip(MGBimage{i},1));
    colormap(gray)
    % xticks([]);
    % yticks([]);
    axis xy
    daspect([1 1 1]);
    hold on

    xlimits = get(gca,'xlim');
    ylimits = get(gca,'ylim');

    I = flip(MGBimage{i},1);

    [C,h]=imcontour(I,1);
    [~, col] = find(C == h.LevelList);
    for j = 1:length(col)
        if j == length(col)
            inner{i}{j} = C(:,(col(j)+1):end);
        else
            inner{i}{j} = C(:,(col(j)+1):col(j+1)-1);
        end
        plot(inner{i}{j}(1,:),inner{i}{j}(2,:),'Color',RegionColor{j});
        
        inner_adj{i}{j}(1,:) = inner{i}{j}(1,:)*(xrange(2)-xrange(1))/size(MGBimage{i},2)+xrange(1);
        inner_adj{i}{j}(2,:) = inner{i}{j}(2,:)*(yrange(2)-yrange(1))/size(MGBimage{i},1)+yrange(1);
    end

end

% Imcontour generates an overall outline as 1st contour, and then draw inner segments in the order of sizes.
% Here, 1st one is excluded, and contours corresponding to MGv, MGd, SG, MGm, BIC, and MZMG are determined manually.
MGBoutlines = cell(6,length(MGBimage));
for i = 1:length(inner_adj)
    if ismember(i,[1,5,6])
        MGBoutlines{1,i} = inner_adj{i}{2}; % MGv
        MGBoutlines{2,i} = inner_adj{i}{5}; % MGd
        MGBoutlines{3,i} = inner_adj{i}{3}; % MGm
        MGBoutlines{4,i} = inner_adj{i}{4}; % SG
    elseif ismember(i,[2,3,4])
        MGBoutlines{1,i} = inner_adj{i}{2}; % MGv
        MGBoutlines{2,i} = inner_adj{i}{4}; % MGd
        MGBoutlines{3,i} = inner_adj{i}{3}; % MGm
        MGBoutlines{4,i} = inner_adj{i}{5}; % SG
    elseif ismember(i,[7])
        MGBoutlines{1,i} = inner_adj{i}{2}; % MGv
        MGBoutlines{3,i} = inner_adj{i}{3}; % MGm (to connect with BIC)
        MGBoutlines{5,i} = inner_adj{i}{3}; % BIC
    elseif ismember(i, [8 ])
        MGBoutlines{1,i} = inner_adj{i}{2}; % MGv ( to connect with MZMG)
        MGBoutlines{6,i} = inner_adj{i}{2}; % MZMG
        MGBoutlines{5,i} = inner_adj{i}{3}; % BIC
    elseif ismember(i, [9])
        MGBoutlines{6,i} = inner_adj{i}{2}; % MZMG
        MGBoutlines{5,i} = inner_adj{i}{3}; % BIC
    elseif ismember(i, [10,11])
        MGBoutlines{5,i} = inner_adj{i}{2}; % BIC
    end
        
end


for i = 1:length(inner_adj)
    MGBdrawdata(i).APposition = APpos(i);
    if ismember(i,[1,5,6])
        MGBdrawdata(i).MGv = inner_adj{i}{2};
        MGBdrawdata(i).MGd = inner_adj{i}{5};
        MGBdrawdata(i).MGm = inner_adj{i}{3};
        MGBdrawdata(i).SG = inner_adj{i}{4};
    elseif ismember(i,[2,3,4])
        MGBdrawdata(i).MGv = inner_adj{i}{2};
        MGBdrawdata(i).MGd = inner_adj{i}{4};
        MGBdrawdata(i).MGm = inner_adj{i}{3};
        MGBdrawdata(i).SG = inner_adj{i}{5};
    elseif ismember(i,[7])
        MGBdrawdata(i).MGv = inner_adj{i}{2};
        MGBdrawdata(i).MGm = inner_adj{i}{3}; % to connect MGm and BIC
        MGBdrawdata(i).BIC = inner_adj{i}{3};
    elseif ismember(i,[8])
        MGBdrawdata(i).MGv = inner_adj{i}{2}; % to connect MGv and MZMG
        MGBdrawdata(i).MZMG = inner_adj{i}{2};
        MGBdrawdata(i).BIC = inner_adj{i}{3};
    elseif ismember(i,[9])
        MGBdrawdata(i).MZMG = inner_adj{i}{2};
        MGBdrawdata(i).BIC = inner_adj{i}{3};
    elseif ismember(i,[10,11])
        MGBdrawdata(i).BIC = inner_adj{i}{2};
        
    end
        
end



%%
f2 = figure('position',[50 200 600 600]);
for i = 1:length(inner_adj)
    hold on
    grid on

    for j = 2:length(inner_adj{i})
        plot3(APpos(i)*ones(1,size(inner_adj{i}{j},2)),inner_adj{i}{j}(1,:),inner_adj{i}{j}(2,:),'Color',RegionColor{j});
    end
end

% daspect([0.001 1 1]);

xlim([-4.5 -2.5]);
set(gca,'Xdir','reverse')
xlabel('Rostral-Caudal (mm, from Bregma)');

yticks(1.5:0.25:2.5);
ylabel('Medial-Lateral (mm)');

zticks(2.5:0.25:4);
set(gca,'Zdir','reverse')
zlabel('Dorsal-Ventral (mm)');
set(gca,'Xcolor',[0 0 0],'Ycolor',[0 0 0],'TickLength',[0.045 0.025]);
view(130,20)


save([filepath filesep 'MGBdrawdata.mat'], 'MGBdrawdata');
