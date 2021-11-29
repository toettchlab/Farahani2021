% Plotting KTR data from SteadyState_Analysis.m
% Payam Farahani

clear all
close all
clc

% Directory containing spreadsheets

addpath('./KTR measurements');

currDir = pwd; 
KTRDir = '/KTR measurements';

load analysis

%% Plot results

% plot temporal heatmaps of C/N ratios
figure(1)
clims = [0 2.5];
imagesc(experiment(1).CN(1:20,1:n_quant),clims)
colormap jet
xlabel('time (min)')
title('C/N ratios')

% plot individual KTR trajectories with Erk pulses labeled
figure(2)
for k = 1:length(imlist)
    for j = 1:length(experiment(k).CN(:,1))
        plot(time,experiment(k).CN(j,:),(experiment(k).peakLoc_matrix(j,1:experiment(k).n_pulse(j))-1)*3, experiment(k).peakMag_matrix(j,1:experiment(k).n_pulse(j)),'o')
        ylim([0 3])
        pause
    end
end
xlabel('time (min)')
ylabel('Erk activity (C/N ratio)')

