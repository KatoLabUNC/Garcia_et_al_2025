% Determine layer identity of individual channels for linear probe recordings in the auditory cortex.
% HK 240406.
%
% Prerequisite:
% 1. determine surface and white matter channels based on spike distributions.
% 2. estimate L4 channels based on current source density analysis.
% 3. generate LFP data. (We downsampled LFP to 250Hz.)
%
% Cross-channel correlation is used to determine layer boundaries.
% This method typically segregates channels into four groups: layer 1, layer 2-4, layer 5, and layer 6. 
% ref: Senzai et al., Neuron 2019.  https://pubmed.ncbi.nlm.nih.gov/30635232/
% We determine the L3-L4 border based on the L4-top channels in CSD signals.
% We determine the L4-L5 border as the mean of CSD- and LFP-based boundaries.
%
% Usually, correlation based on broadband LFP data works well, but in some experiments, gamma-band or delta-band LFP gave better layer segregation. 
% This code shows cross-channel correlation for all different bands.
%
%

clear all

parentdir = '(parent directory path)';

mouse = '(mouse)';
date = '(date)';
site = '(site ID)';

areaname = 'A2'; % A1 or A2

% pre-determined data
surfaceCh = 10;
whitematterCh = 60;
L4Ch_CSD = [27 31]; % In this example, we used 64-channel linear probe recording, with channel 1 outside the brain and channel 64 within the white matter.


subdir = [parentdir filesep date filesep mouse filesep site];

% load your own LFP data. This can be either spontaneous activity measurement or sound presentation experiment.
load('(your LFP file path)');


%% store data based on clicks CSD

layer_info = [];
layer_info.mouse = mouse;
layer_info.date = date;
layer_info.site = site;
layer_info.areaname = areaname;
layer_info.surfaceCh = surfaceCh;
layer_info.whitematterCh = whitematterCh;
layer_info.L4Ch_CSD = L4Ch_CSD;

            
%% calculate across-channels correlation

% just use channels within the cortex

curr_data = downsample_traces;
curr_data([1:surfaceCh, whitematterCh:end],:) = [];
L4_adj = L4Ch_CSD-surfaceCh;
ch_num = whitematterCh-surfaceCh-1;


h = figure('position',[100 100 1250 285]);

subplot1(1,5,'Min',[0.1 0.1],'Max',[0.95 0.8],'Gap',[0.02 0.02]);

climits = [0.3 1];


% plot correlation coefficient across channels.
[b1 a1] = butter(3, 64/downsampled_samplerate*2, 'low');
[b2 a2] = butter(3, 1/downsampled_samplerate*2, 'high');
data_broadband = (filtfilt(b1,a1,curr_data'))';
data_broadband = (filtfilt(b2,a2,data_broadband'))';

R = corrcoef(data_broadband');
subplot1(1); hold on;
imagesc(R);
line([L4_adj(1) L4_adj(1)],[0 ch_num],'color','k');
line([0 ch_num],[L4_adj(1) L4_adj(1)],'color','k');
line([L4_adj(2) L4_adj(2)],[0 ch_num],'color','k');
line([0 ch_num],[L4_adj(2) L4_adj(2)],'color','k');
    xlabel('Channel numb (pia to wm)');
    title('Full trace');
ylabel({[mouse ' Surf ' num2str(surfaceCh)], 'Channel numb (pia to wm)'});
xlim([0 ch_num]);
ylim([0 ch_num]);
set(gca,'ydir','reverse','XaxisLocation','top');
set(gca,'clim',climits);
            
            
% correlation coefficient for traces with different frequency bands

% % delta oscillation

[b1 a1] = butter(3, 4/downsampled_samplerate*2, 'low');
[b2 a2] = butter(3, 1/downsampled_samplerate*2, 'high');
data_delta = (filtfilt(b1,a1,curr_data'))';
data_delta = (filtfilt(b2,a2,data_delta'))';

R_delta = corrcoef(data_delta');
subplot1(2); hold on;
imagesc(R_delta);
line([L4_adj(1) L4_adj(1)],[0 ch_num],'color','k');
line([0 ch_num],[L4_adj(1) L4_adj(1)],'color','k');
line([L4_adj(2) L4_adj(2)],[0 ch_num],'color','k');
line([0 ch_num],[L4_adj(2) L4_adj(2)],'color','k');
    xlabel('Channel numb (pia to wm)');
    title('Delta band');
xlim([0 ch_num]);
ylim([0 ch_num]);
set(gca,'ydir','reverse','XaxisLocation','top');
set(gca,'clim',climits);
            
% % theta to alpha oscillation

[b1 a1] = butter(3, 12/downsampled_samplerate*2, 'low');
[b2 a2] = butter(3, 4/downsampled_samplerate*2, 'high');
data_theta = (filtfilt(b1,a1,curr_data'))';
data_theta = (filtfilt(b2,a2,data_theta'))';


R_theta = corrcoef(data_theta');
subplot1(3); hold on;
imagesc(R_theta);
line([L4_adj(1) L4_adj(1)],[0 ch_num],'color','k');
line([0 ch_num],[L4_adj(1) L4_adj(1)],'color','k');
line([L4_adj(2) L4_adj(2)],[0 ch_num],'color','k');
line([0 ch_num],[L4_adj(2) L4_adj(2)],'color','k');
    xlabel('Channel numb (pia to wm)');
    title('Theta-Alpha band');
set(gca,'ydir','reverse','XaxisLocation','top');
set(gca,'clim',climits);
xlim([0 ch_num]);
ylim([0 ch_num]);
            
% % beta to gamma oscillation

[b1 a1] = butter(3, 64/downsampled_samplerate*2, 'low');
[b2 a2] = butter(3, 12/downsampled_samplerate*2, 'high');
data_gamma = (filtfilt(b1,a1,curr_data'))';
data_gamma = (filtfilt(b2,a2,data_gamma'))';

R_gamma = corrcoef(data_gamma');
subplot1(4); hold on;
imagesc(R_gamma);
line([L4_adj(1) L4_adj(1)],[0 ch_num],'color','k');
line([0 ch_num],[L4_adj(1) L4_adj(1)],'color','k');
line([L4_adj(2) L4_adj(2)],[0 ch_num],'color','k');
line([0 ch_num],[L4_adj(2) L4_adj(2)],'color','k');
    xlabel('Channel numb (pia to wm)');
    title('Beta-Gamma band');
set(gca,'ydir','reverse','XaxisLocation','top');
set(gca,'clim',climits);
xlim([0 ch_num]);
ylim([0 ch_num]);
            

% try k-means clustering on broadband
initialcenters = [mean(data_broadband(1:round(L4_adj(1)/2),:),1);... % upper half of supragranular
    mean(data_broadband(L4_adj(1):L4_adj(2),:),1);... % CSD L4
    mean(data_broadband(L4_adj(2)+1:round((L4_adj(2)+ch_num)/2),:),1);... % upper half of subgranular
    mean(data_broadband(round((L4_adj(2)+ch_num)/2)+1:ch_num,:),1)];
idx = kmeans(data_broadband,4,'Distance','correlation','start',initialcenters);

subplot1(5); hold on;
plot(idx,1:ch_num,'k');
set(gca,'ydir','reverse','box','off');
xlim([0.975 4.025]);
ylim([0 ch_num]);
    title('k-means clustering');
        

% % try k-means clustering on gamma band
% initialcenters = [mean(data_gamma(1:round(L4_adj(1)/2),:),1);... % upper half of supragranular
%     mean(data_gamma(L4_adj(1):L4_adj(2),:),1);... % CSD L4
%     mean(data_gamma(L4_adj(2)+1:round((L4_adj(2)+ch_num)/2),:),1);... % upper half of subgranular
%     mean(data_gamma(round((L4_adj(2)+ch_num)/2)+1:ch_num,:),1)];
% idx = kmeans(data_gamma,4,'Distance','correlation','start',initialcenters);
% 
% subplot1(5); hold on;
% plot(idx,1:ch_num,'k');
% set(gca,'ydir','reverse','box','off');
% xlim([0.975 4.025]);
% ylim([0 ch_num]);
% title('k-means clustering');
        
% % try k-means clustering on delta band
% initialcenters = [mean(data_delta(1:round(L4_adj(1)/2),:),1);... % upper half of supragranular
%     mean(data_delta(L4_adj(1):L4_adj(2),:),1);... % CSD L4
%     mean(data_delta(L4_adj(2)+1:round((L4_adj(2)+ch_num)/2),:),1);... % upper half of subgranular
%     mean(data_delta(round((L4_adj(2)+ch_num)/2)+1:ch_num,:),1)];
% idx = kmeans(data_delta,4,'Distance','correlation','start',initialcenters);
% 
% subplot1(5); hold on;
% plot(idx,1:ch_num,'k');
% set(gca,'ydir','reverse','box','off');
% xlim([0.975 4.025]);
% ylim([0 ch_num]);
% title('k-means clustering');
        

idx_smo = round(smooth(idx,5));
plot(idx_smo,1:ch_num,'r');
        

L2_top = surfaceCh + find(idx_smo==2,1,'first'); % note, this likely includes L2/3. Not used.
L4_bottom = surfaceCh + find(idx_smo==2,1,'last');
L5_bottom = surfaceCh + find(idx_smo==3,1,'last');
        
text(2.5,8,['surface: ' num2str(surfaceCh)],'fontsize',10);      
text(2.5,18,['L2/3 top: ' num2str(L2_top)],'fontsize',10);      
text(2.5,28,['L4 bottom: ' num2str(L4_bottom)],'fontsize',10);      
text(2.5,38,['L5 bottom: ' num2str(L5_bottom)],'fontsize',10);      


%% check and edit as necessary

% overwrite the layer_info file?
approved = questdlg('Layer borders OK?');

if strcmp(approved, 'No')
    f = msgbox('Manually edit L2_top, L4_bottom, and L5_bottom. You may try using different frequency bands for k-means clustering. After manual edits, continue with the rest of the code.');
    error('Please edit area borders.');
end

% check the visualized plot and edit as needed.
% k-means clustering usually works well, but sometimes, L1-L2 boundary may generate noise.

layer_info.L2top_LFP = L2_top;
layer_info.L4bottom_LFP = L4_bottom;
layer_info.L5bottom_LFP = L5_bottom;
layer_info.totalCh = size(downsample_traces,1);



new_L4_bottom = round(mean([layer_info.L4Ch_CSD(2), layer_info.L4bottom_LFP]));

layer_info.ch_label = [];

% layer 1
for ch = layer_info.surfaceCh+1:layer_info.L2top_LFP-1
    layer_info.ch_label{ch} = [layer_info.areaname 'L1'];
end

% layer 2/3
for ch = layer_info.L2top_LFP:layer_info.L4Ch_CSD(1)-1
    layer_info.ch_label{ch} = [layer_info.areaname 'L23'];
end

% layer 4
for ch = layer_info.L4Ch_CSD(1):new_L4_bottom
    layer_info.ch_label{ch} = [layer_info.areaname 'L4'];
end

% layer 5
for ch = new_L4_bottom+1:layer_info.L5bottom_LFP
    layer_info.ch_label{ch} = [layer_info.areaname 'L5'];
end

% layer 6
for ch = layer_info.L5bottom_LFP+1:layer_info.whitematterCh-1
    layer_info.ch_label{ch} = [layer_info.areaname 'L6'];
end

% below cortex
for ch = layer_info.whitematterCh:layer_info.totalCh
    layer_info.ch_label{ch} = [];
end

save([subdir filesep 'layer_info.mat'],'layer_info');





