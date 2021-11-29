% Analysis of KTR using C/N ratio
% Payam Farahani

clear all
close all
clc

addpath('./utils');

%% Enter parameters

n_quant = 121; % time over which Erk activity is averaged
time_step = 3; % min between each acquisition;
time = [0:time_step:(n_quant-1)*time_step]; % time vector

%% Directory settings

% Directory containing spreadsheets
currDir = pwd; 
KTRDir = '/KTR measurements';

% Import data and list of tissues
cd(strcat(currDir,KTRDir));
dirNAME = dir;
namArr = {'2019_02_07_PF2_130Pa_well1.xlsx'; 
        '2019_02_07_PF2_hydrogel_130Pa_well2.xlsx'; 
        '2020_01_08_PF1_130Pa'};
imlist = namArr;

% Generate data structure containing experiment information
for j = 1:length(imlist)
    explist = char(imlist(j));
    experiment(j).pos = char(explist(3:(end-5)));
    experiment(j).file = explist;
end

% clear explist imlist namArr dirNAME
cd(currDir);

%% Read time values and C/N ratios for each experiment

for j = 1:size(experiment,2)
    cd(strcat(currDir,KTRDir));
    
    % read xlsx file
    raw_vals = xlsread(experiment(j).file,1,'A4:CB124'); % raw C and N values
    experiment(j).raw_vals = raw_vals(:,all(~isnan(raw_vals)));
    
    experiment(j).BG = xlsread(experiment(j).file,5,'F2'); % background intensity
    
   % subtract background from raw intensities
    experiment(j).vals_noBG = experiment(j).raw_vals - experiment(j).BG;
    
   % calculate C/N ratios
   experiment(j).CN = [];
   
   for k = 1:length(experiment(j).raw_vals(1,:))/2
       tempCN = experiment(j).vals_noBG(:,2*k-1) ./ experiment(j).vals_noBG(:,2*k);
       experiment(j).CN = [experiment(j).CN tempCN];
   end
  
   experiment(j).CN = experiment(j).CN';
   
   % running average of C/N to eliminate internal noise
   for k = 1:length(experiment(j).CN(:,1))
       experiment(j).CN(k,1) = mean(experiment(j).CN(k,1:2));
       for l = 2:length(experiment(j).CN(k,:))-1
           experiment(j).CN(k,l) = mean(experiment(j).CN(k,l-1:l+1));
       end
       experiment(j).CN(k,end) = mean(experiment(j).CN(k,end-1:end));
   end
   
   nCells = length(experiment(j).CN(:,1));
   
   % calculate time-averaged C/N for each cell
   for k = 1:nCells
       experiment(j).CN_time_avg(k) = mean(experiment(j).CN(k,:));
   end
   
   % quantify pulses
   sel = 0.2; % the amount above surrounding data for a peak to be
   thresh = 0.5; % threshold value which peaks must be larger than to be maxima
   extrema = 1; % 1 if maxima are desired, -1 if minima are desired
   includeEndpoints = false; % If true the endpoints will be included as possible extrema otherwise they will not be included
   interpolate = false;
    
   peakLoc_matrix = zeros(nCells,100);
   peakMag_matrix = zeros(nCells,100);
    
   for k = 1:nCells
       [peakLoc, peakMag] = peakfinder(experiment(j).CN(k,1:n_quant), sel, thresh, extrema, includeEndpoints, interpolate);
       experiment(j).n_pulse(k) = length(peakLoc);
       experiment(j).peakLoc_matrix(k,1:length(peakLoc)) = peakLoc;
       experiment(j).peakMag_matrix(k,1:length(peakMag)) = peakMag;
   end
    
    experiment(j).n_pulse = experiment(j).n_pulse';
  
    
   % quantify proportion of cells in each dynamic state
   experiment(j).dyn_num = zeros(1,3);
   for k = 1:nCells
       if experiment(j).n_pulse(k) < 2
           if experiment(j).CN_time_avg(k) < 1
               experiment(j).dyn_num(1) = experiment(j).dyn_num(1) + 1; % constantly off
           else
               experiment(j).dyn_num(3) = experiment(j).dyn_num(3) + 1; % constantly on
           end
       elseif experiment(j).n_pulse(k) > 1
           experiment(j).dyn_num(2) = experiment(j).dyn_num(2) + 1;
       end
   end

    experiment(j).dyn_frac = experiment(j).dyn_num / nCells;


    % quantify proportion of cells with x number of pulses
    
    [experiment(j).pulse_histogram(:,1), experiment(j).pulse_histogram(:,2)] = groupcounts(experiment(j).n_pulse);
    experiment(j).pulse_histogram(:,1) = experiment(j).pulse_histogram(:,1) ./ length(experiment(j).n_pulse);

end

% compile time-average C/N and pulses for all experiments
all_CN_time_avg = [];
all_n_pulse = [];

for j = 1:length(imlist)
    all_CN_time_avg = [all_CN_time_avg; experiment(j).CN_time_avg'];
    all_n_pulse = [all_n_pulse; experiment(j).n_pulse];
end

% calculate mean and std dev for dynamic demographics
all_dyn_frac = [];

for j = 1:size(imlist)
    all_dyn_frac = [all_dyn_frac; experiment(j).dyn_frac];
end

mean_dyn_frac = mean(all_dyn_frac,1);
std_dyn_frac = std(all_dyn_frac,1);

% calculate mean and std dev for fraction of cells with x pulses

pulse_histogram = zeros(10,size(experiment,2));

for j = 1:size(experiment,2)
    pulse_histogram(1:length(experiment(j).pulse_histogram(:,2)),j) = experiment(j).pulse_histogram(:,1);
end

pulse_histogram_mean = zeros(10,1);
pulse_histogram_std = zeros(10,1);
for j = 1:length(pulse_histogram(:,1))
    pulse_histogram_mean(j) = mean(pulse_histogram(j,:));
    pulse_histogram_std(j) = std(pulse_histogram(j,:));
end

% calculate population-averaged C/N for every time point
all_CN = [];

for j = 1:size(experiment,2)
    all_CN = [all_CN; experiment(j).CN];
end

for k = 1:length(time)
    CN_pop_avg(k,1) = mean(all_CN(:,k)); % population-averaged C/N
    CN_pop_avg(k,2) = std(all_CN(:,k)); % standard deviation of C/N
end


%% Save results

filename = 'analysis.mat';
save(filename)