% Payam Farahani
% Script for analyzing puncta from z-stack images

clear all
close all
clc

%% Enter parameters

% EGF object size thresholds (pixels)
minsize_EGF = 10;
maxsize_EGF = 200;

% EGFR object size thresholds (pixels)
minsize_EGFR = 10;
maxsize_EGFR = 200;

% image scaling
dim_xy = 1/4.2221; % µm/pixel in xy plane
dim_z = 0.3; % µm z-step

% sigma values for Gaussian filters
sigma_EGF = 5;
sigma_EGFR = 5;

%% Get file info

% get file names
files = dir('*_EGF.tif');

for n = 1:length(files)
    Data(n).file = files(n).name(1:end-8);

% name images
name_EGF = append(Data(n).file,'_EGF.tif');
name_EGFR = append(Data(n).file,'_EGFR.tif');

% get file info
info_EGF = imfinfo(name_EGF);
info_EGFR = imfinfo(name_EGFR);

n_slices = length(info_EGF); % number of slices in z-stack

% initialize image arrays: import raw images
EGF = zeros(info_EGF(1).Width(1),info_EGF(1).Width(1),n_slices);
EGFR = EGF;

for j = 1:n_slices
    EGF(:,:,j) = im2double(imread(name_EGF,j));
    EGFR(:,:,j) = im2double(imread(name_EGFR,j));
end

%% Background subtraction

% apply wide Gaussian filter to image then subtract from raw image

% initialize filtered image arrays
% fnuclei = nuclei.*0;
fEGF = EGF.*0;
fEGFR = fEGF;

% subtract Gaussian blurs from raw images
for j = 1:n_slices
    gauss_EGF = imgaussfilt(EGF(:,:,j),sigma_EGF); 
    fEGF(:,:,j) = EGF(:,:,j) - gauss_EGF;
    
    gauss_EGFR = imgaussfilt(EGFR(:,:,j),sigma_EGFR);
    fEGFR(:,:,j) = EGFR(:,:,j) - gauss_EGFR;
end

%% Set intensity threshold for EGF puncta

threshold_EGF = mean(nonzeros(fEGF(:))) + 6*std(nonzeros(fEGF(:)));
threshold_EGFR = mean(nonzeros(fEGFR(:))) + 6*std(nonzeros(fEGFR(:)));

Data(n).treshold_EGF = threshold_EGF;
Data(n).treshold_EGFR = threshold_EGFR;

%% Threshold and count particles

% select pixels greater than intensity threshold
fEGF = fEGF > threshold_EGF;
fEGFR = fEGFR > threshold_EGFR;

% select pixels positive for both EGF and EGFR
Coloc = fEGF.*fEGFR;

% filter objects below minimum size
fEGF = bwareaopen(fEGF,minsize_EGF);
fEGFR = bwareaopen(fEGFR,minsize_EGFR);
Coloc = bwareaopen(Coloc,minsize_EGF);

% filter objects above maximum size
CC_EGF = bwconncomp(fEGF);
S_EGF = regionprops3(CC_EGF, 'Volume');
L_EGF = labelmatrix(CC_EGF);
fEGF = ismember(L_EGF, find([S_EGF.Volume] <= maxsize_EGF));

CC_EGFR = bwconncomp(fEGFR);
S_EGFR = regionprops3(CC_EGFR, 'Volume');
L_EGFR = labelmatrix(CC_EGFR);
fEGFR = ismember(L_EGFR, find([S_EGFR.Volume] <= maxsize_EGFR));

CC_Coloc = bwconncomp(Coloc);
S_Coloc = regionprops3(CC_Coloc, 'Volume');
L_Coloc = labelmatrix(CC_Coloc);
Coloc = ismember(L_Coloc, find([S_Coloc.Volume] <= maxsize_EGF));

% identify objects
cc_EGF = bwconncomp(fEGF);
r_EGF = regionprops(cc_EGF,'PixelList');
Data(n).volumes_EGF = zeros(cc_EGF.NumObjects,1);

cc_EGFR = bwconncomp(fEGFR);
r_EGFR = regionprops(cc_EGFR,'PixelList');
Data(n).volumes_EGFR = zeros(cc_EGFR.NumObjects,1);

cc_Coloc = bwconncomp(Coloc);
r_Coloc = regionprops(cc_Coloc,'PixelList');
Data(n).volumes_Coloc = zeros(cc_Coloc.NumObjects,1);


% calculate volumes of objects
for k = 1:cc_EGF.NumObjects
    Data(n).volumes_EGF(k) = length(r_EGF(k).PixelList)*dim_xy^2*dim_z; % compute volumes from pixel number
end

for k = 1:cc_EGFR.NumObjects
    Data(n).volumes_EGFR(k) = length(r_EGFR(k).PixelList)*dim_xy^2*dim_z; % compute volumes from pixel number
end

for k = 1:cc_Coloc.NumObjects
    Data(n).volumes_Coloc(k) = length(r_Coloc(k).PixelList)*dim_xy^2*dim_z; % compute volumes from pixel number
end

Data(n).volume_total_EGF = sum(Data(n).volumes_EGF);
Data(n).volume_mean_EGF = mean(Data(n).volumes_EGF);
Data(n).number_of_particles_EGF = cc_EGF.NumObjects;

Data(n).volume_total_EGFR = sum(Data(n).volumes_EGFR);
Data(n).volume_mean_EGFR = mean(Data(n).volumes_EGFR);
Data(n).number_of_particles_EGFR = cc_EGFR.NumObjects;

Data(n).volume_total_Coloc = sum(Data(n).volumes_Coloc);
Data(n).volume_mean_Coloc = mean(Data(n).volumes_Coloc);
Data(n).number_of_particles_Coloc = cc_Coloc.NumObjects;    
        
% generate test files
test_file_name_EGF = sprintf([name_EGF(1:end-4) '_seg.tif']);
for k = 1:n_slices
    imwrite(fEGF(:,:,k),test_file_name_EGF,'WriteMode','append','compression','none');
end

test_file_name_EGFR = sprintf([name_EGFR(1:end-4) '_seg.tif']);
for k = 1:n_slices
    imwrite(fEGFR(:,:,k),test_file_name_EGFR,'WriteMode','append','compression','none');
end

test_file_name_Coloc = sprintf([name_EGF(1:end-4) '_coloc.tif']);
for k = 1:n_slices
    imwrite(Coloc(:,:,k),test_file_name_Coloc,'WriteMode','append','compression','none');
end

end

% save .mat file
filename = 'Data.mat';
save(filename)



