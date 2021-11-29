% Payam Farahani

clear all
close all
clc

%% Enter parameters

minsize_EGF = 2; % minimum pixel size for EGF puncta
maxsize_EGF = 200; % maximum pixel size for EGF puncta

dim_xy = 1/9.2308; % scaling in µm/pixel

EGFR_BG = 120;

%% Get file info

% get file names
name_EGF = uigetfile('*.tif','EGF'); % max projection of EGF-488 channel (rolling ball BG subtracted + gaussian filtered)
name_EGFR = uigetfile('*.tif','EGFR'); % max projection of EGFR channel with no BG subtraction

% get file info
info_EGF = imfinfo(name_EGF);
info_EGFR = imfinfo(name_EGFR);

% initialize image arrays
EGF = zeros(info_EGF(1).Height(1),info_EGF(1).Width(1));
EGF = im2double(imread(name_EGF));
EGFR = zeros(info_EGFR(1).Height(1),info_EGFR(1).Width(1));
EGFR = double(imread(name_EGFR)) - EGFR_BG; %subtract BG fluorescence

%% Set intensity threshold for EGF puncta

threshold_EGF = mean(EGF(:)) + 1*std(EGF(:));

%% Threshold EGF-488 images

% filter image by intensity threshold
fEGF = EGF > threshold_EGF;
fEGF_intfilt = fEGF;

% filter objects below minimum size
fEGF = bwareaopen(fEGF,minsize_EGF);
fEGF_minfilt = fEGF;

% filter objects above maximum size
CC = bwconncomp(fEGF);
S = regionprops(CC, 'Area');
L = labelmatrix(CC);
fEGF = ismember(L, find([S.Area] <= maxsize_EGF));
fEGF_maxfilt = fEGF;

% identify objects
cc_fEGF = bwconncomp(fEGF);
r_fEGF = regionprops(cc_fEGF,'PixelList');

%% Create masks of cells

run EGF_Membrane_Outlines

%% Assign puncta to cells

% coordinates of puncta in each cell

for j = 1:length(data)
    data(j).puncta = [];
    data(j).unfilt = [];
end

for j = 1:length(r_fEGF);
    for k = 1:length(data)
        test = inpolygon(r_fEGF(j).PixelList(:,1),r_fEGF(j).PixelList(:,2),data(k).coords(:,1),data(k).coords(:,2));
        data(k).puncta = [data(k).puncta; r_fEGF(j).PixelList(test==1,:)];
    end
end

%% Unfiltered signal per cell

% identify objects
cc_EGF = bwconncomp(EGF);
r_EGF = regionprops(cc_EGF,'PixelList');

for j = 1:length(r_EGF)
    for k = 1:length(data)
        test = inpolygon(r_EGF(j).PixelList(:,1),r_EGF(j).PixelList(:,2),data(k).coords(:,1),data(k).coords(:,2));
        data(k).unfilt = [data(k).unfilt; r_EGF(j).PixelList(test==1,:)];
    end
end

%% Analyze puncta

for j = 1:length(data)
    
    % area of puncta per cell
    data(j).areapuncta = length(data(j).puncta) * dim_xy^2;
    
    % puncta area / cell area (per cell)
    data(j).areapuncta_per_areacell = data(j).areapuncta / data(j).areas;
    
    % unfiltered intensity per cell
    for k = 1:length(data(j).unfilt)
        data(j).unfilt_int(k) = EGF(data(j).unfilt(k,2), data(j).unfilt(k,1));
    end
    
    % sum of unfiltered pixels per cell
    data(j).unfilt_sum = sum(data(j).unfilt_int);
       
end

%% Analyze EGFR

for j = 1:length(data)
    data(j).EGFR_pixels = double(poly2mask(data(j).coords(:,1),data(j).coords(:,2),rows,cols)) .* EGFR;
    data(j).EGFR_sum = sum(data(j).EGFR_pixels(:));
    data(j).EGFR_percell = data(j).EGFR_sum ./ data(j).areas;
end

%% Output and saving

% generate test files
test_file_name_EGF = sprintf([name_EGF(1:end-4) '_seg.tif']);
imwrite(fEGF(:,:),test_file_name_EGF,'WriteMode','overwrite','compression','none')

% save .mat file
filename = sprintf(name_EGF(1:end-8),'.mat');
save(filename)

%% Validation

figure(1)
subplot(1,3,1)
imshow(EGF)
subplot(1,3,2)
imshow(fEGF)
subplot(1,3,3)
imshow(masks)

figure(2)
testim = zeros(rows,cols);
for l = 1:length(data)
    tempmask = data(l).coords;
    temppuncta = data(l).puncta;
    
    for j = 1:length(tempmask)
        testim(tempmask(j,2),tempmask(j,1)) = l;
    end
    
    for k = 1:length(temppuncta)
        testim(temppuncta(k,2),temppuncta(k,1)) = l;
    end
end
imshow(testim)
caxis([0 length(data)])
colormap colorcube

