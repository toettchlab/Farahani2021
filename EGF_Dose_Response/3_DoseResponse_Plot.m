% Plotting KTR data from DoseResponse_Analysis.m
% Payam Farahani

clear all
close all
clc

% Directory containing spreadsheets

addpath('./KTR measurements');

currDir = pwd; 
KTRDir = '/KTR measurements';

load analysis

%% Plot temporal heatmaps of C/N ratios

figure(1)

for j = 1:length(pooled)
    subplot(3,3,j)
    clims = [0.3 2.5];
    imagesc(pooled(j).CN_sort,clims)
end

colormap jet
xlabel('time (min)')
title('C/N ratios')

%% Plot population-averaged C/N vs. time

figure(2)
subplot(1,3,1)
hold on
for j = 1:3
    plot(pooled(j).time,pooled(j).CN_pop_avg,'-')
end
ylim([0 1.8])
hold off

subplot(1,3,2)
hold on
for j = 4:6
    plot(pooled(j).time,pooled(j).CN_pop_avg,'-')
end
ylim([0 1.8])
hold off

subplot(1,3,3)
hold on
for j = 7:9
    plot(pooled(j).time,pooled(j).CN_pop_avg,'-')
end
ylim([0 1.8])
hold off
legend('130 Pa','910 Pa','4020 Pa')

%% Plot area-under-curve of C/N trajectories

figure(3)

subplot(1,2,1)
hold on
for j = 1:3:7
    plot(pooled(j).CN_int,pooled(j).CN_peak,'o')
    plot(pooled(j+1).CN_int,pooled(j+1).CN_peak,'s')
    plot(pooled(j+2).CN_int,pooled(j+2).CN_peak,'d')
end
hold off
xlabel('C/N AUC')
ylabel('peak C/N')

subplot(1,2,2)
hold on
for j = 1:3:7
    plot(pooled(j).mean_CN_end,pooled(j).mean_CN_peak,'o')
    plot(pooled(j+1).mean_CN_end,pooled(j+1).mean_CN_peak,'s')
    plot(pooled(j+2).mean_CN_end,pooled(j+2).mean_CN_peak,'d')
end
hold off
legend('130 Pa (0.2 ng/mL)','910 Pa (0.2 ng/mL)','4020 Pa (0.2 ng/mL)','130 Pa (2 ng/mL)','910 Pa (2 ng/mL)','4020 Pa (2 ng/mL)','130 Pa (20 ng/mL)','910 Pa (20 ng/mL)','4020 Pa (20 ng/mL)')
xlabel('final C/N')
ylabel('peak C/N')
xlim([0 1])
ylim([0 2])

%% Box plots

CN_peak = [];
CN_int = [];
CN_end = [];
group = [];

for j = 1:length(pooled)
    CN_peak = [CN_peak pooled(j).CN_peak'];
    CN_int = [CN_int pooled(j).CN_int'];
    CN_end = [CN_end pooled(j).CN_end'];
    group = [group ones(1,length(pooled(j).CN_peak))*j];
end

figure(5)
subplot(1,3,1)
boxplot(CN_peak,group)
title('peak C/N ratio')
subplot(1,3,2)
boxplot(CN_int,group)
title('AUC C/N ratio')
subplot(1,3,3)
boxplot(CN_end,group)
title('final C/N ratio')

%% heatmaps

% peak C/N
for j = 1:length(pooled)
    heatmap_CN_peak(j) = pooled(j).mean_CN_peak;
    heatmap_CN_end(j) = pooled(j).mean_CN_end;
    heatmap_CN_int(j) = pooled(j).mean_CN_int;
end

heatmap_CN_peak = reshape(heatmap_CN_peak,[3,3]);
heatmap_CN_end = reshape(heatmap_CN_end,[3,3]);
heatmap_CN_int = reshape(heatmap_CN_int,[3,3]);

figure(6)
colormap summer
subplot(1,3,1)
imagesc(heatmap_CN_peak)
subplot(1,3,2)
imagesc(heatmap_CN_end)
subplot(1,3,3)
imagesc(heatmap_CN_int)


