% C/N ratio analysis for EGF dose response experiments
% Payam Farahani


clear all
close all
clc

addpath('./utils');

%% Enter parameters

timestep = 1; % minutes between each acquisition

%% Directory settings

% Directory containing spreadsheets
currDir = pwd; 
KTRDir = '/KTR measurements';

% Import data and list of tissues
cd(strcat(currDir,KTRDir));
dirNAME = dir;
namArr = {dirNAME.name}';
imlist = namArr(~cellfun('isempty', regexp(namArr,'xlsx')));

% Generate data structure containing experiment information
for j = 1:length(imlist)
    explist = char(imlist(j));
    experiment(j).file = explist;
    experiment(j).EGF = char(explist(1:3));
    experiment(j).stiffness = char(explist(5:8));
end

% clear explist imlist namArr dirNAME
cd(currDir);

%% Read time values and C/N ratios for each experiment

for j = 1:size(experiment,2)
    cd(strcat(currDir,KTRDir));
    
    % read xlsx file
    raw_vals = xlsread(experiment(j).file,1,'A4:CB254'); % raw C and N values
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
   
   % create time vector
   experiment(j).time = [0:timestep:timestep*(length(experiment(j).CN(1,:))-1)];
   
   % running average of C/N to eliminate internal noise
   for k = 1:length(experiment(j).CN(:,1))
       experiment(j).CN(k,1) = mean(experiment(j).CN(k,1:2));
       for l = 2:length(experiment(j).CN(k,:))-1
           experiment(j).CN(k,l) = mean(experiment(j).CN(k,l-1:l+1));
       end
       experiment(j).CN(k,end) = mean(experiment(j).CN(k,end-1:end));
   end
   
   experiment(j).nCells = length(experiment(j).CN(:,1));
   
   % calculate population-averaged C/N at each time point
   experiment(j).CN_pop_avg = mean(experiment(j).CN);
   
   % calculate time-averaged C/N for each cell
   for k = 1:experiment(j).nCells
       experiment(j).CN_time_avg(k) = mean(experiment(j).CN(k,:));
   end
   
end
   
%% Combine results from each replicate

pooled(1).EGF = '002';
pooled(2).EGF = '002';
pooled(3).EGF = '002';
pooled(4).EGF = '020';
pooled(5).EGF = '020';
pooled(6).EGF = '020';
pooled(7).EGF = '200';
pooled(8).EGF = '200';
pooled(9).EGF = '200';

pooled(1).stiffness = '0130';
pooled(4).stiffness = '0130';
pooled(7).stiffness = '0130';
pooled(2).stiffness = '0910';
pooled(5).stiffness = '0910';
pooled(8).stiffness = '0910';
pooled(3).stiffness = '4020';
pooled(6).stiffness = '4020';
pooled(9).stiffness = '4020';

for j = length(pooled)
    pooled(j).CN = [];
end

for j = 1:length(experiment)
    if contains(experiment(j).EGF,'002') == 1 && contains(experiment(j).stiffness,'0130') == 1
        ind = 1;
    elseif contains(experiment(j).EGF,'002') == 1 && contains(experiment(j).stiffness,'0910') == 1
        ind = 2;
    elseif contains(experiment(j).EGF,'002') == 1 && contains(experiment(j).stiffness,'4020') == 1
        ind = 3;
    elseif contains(experiment(j).EGF,'020') == 1 && contains(experiment(j).stiffness,'0130') == 1
        ind = 4;
    elseif contains(experiment(j).EGF,'020') == 1 && contains(experiment(j).stiffness,'0910') == 1
        ind = 5;
    elseif contains(experiment(j).EGF,'020') == 1 && contains(experiment(j).stiffness,'4020') == 1
        ind = 6;
    elseif contains(experiment(j).EGF,'200') == 1 && contains(experiment(j).stiffness,'0130') == 1
        ind = 7;
    elseif contains(experiment(j).EGF,'200') == 1 && contains(experiment(j).stiffness,'0910') == 1
        ind = 8;
    elseif contains(experiment(j).EGF,'200') == 1 && contains(experiment(j).stiffness,'4020') == 1
        ind = 9;
    end
    
    pooled(ind).CN = [pooled(ind).CN; experiment(j).CN];
end

%% Analyze pooled data

for j = 1:length(pooled)
    
   % create time vector
   pooled(j).time = [0:timestep:timestep*(length(pooled(j).CN(1,:))-1)];
    
   % population-averaged C/N ratio
    pooled(j).CN_pop_avg = mean(pooled(j).CN);
    
   % calculate peak amplitude for each cell
   pooled(j).CN_peak = max(pooled(j).CN(:,1:60),[],2);
   pooled(j).mean_CN_peak = mean(pooled(j).CN_peak);
   pooled(j).std_CN_peak = std(pooled(j).CN_peak);
   
   % calculate final C/N ratio of each cell
   pooled(j).CN_end = mean(pooled(j).CN(:,end-10:end),2);
   pooled(j).mean_CN_end = mean(pooled(j).CN_end);
   pooled(j).std_CN_end = std(pooled(j).CN_end);
   
   % C/N area under curve
   
   % subtract initial C/N from all frames
   pooled(j).CN_corr = pooled(j).CN*0;
   for k = 1:length(pooled(j).CN(:,1))
       pooled(j).CN_corr(k,:) = pooled(j).CN(k,:) - pooled(j).CN(k,1);
   end
   
   % C/N AUC for all time points
   pooled(j).CN_int = sum(pooled(j).CN_corr,2);
   pooled(j).mean_CN_int = mean(pooled(j).CN_int);
   pooled(j).std_CN_int = std(pooled(j).CN_int);
   
   % C/N AUC for first hour of stimulation
   pooled(j).CN_int_early = sum(pooled(j).CN_corr(:,11:70),2);
   pooled(j).mean_CN_int_early = mean(pooled(j).CN_int_early);
   pooled(j).std_CN_int_early = std(pooled(j).CN_int_early);
   
   % C/N AUC after first hour of stimulation
   pooled(j).CN_int_late = sum(pooled(j).CN_corr(:,71:end),2);
   pooled(j).mean_CN_int_late = mean(pooled(j).CN_int_late);
   pooled(j).std_CN_int_late = std(pooled(j).CN_int_late);
   
   % sort CN matrix based on C/N area under curve
   temp = [pooled(j).CN_int pooled(j).CN];
   temp = sortrows(temp,1);
   pooled(j).CN_sort = temp(:,2:end);
   
end


%% Save results

filename = 'analysis.mat';
save(filename)