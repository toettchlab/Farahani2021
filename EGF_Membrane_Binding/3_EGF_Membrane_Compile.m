% Payam Farahani
% combine .mat files for data analysis
%%
clear all
close all
clc

files=dir(strcat('*.mat'));
nf = length(files);

areapuncta.soft = [];
areapuncta.int = [];
areapuncta.stiff = [];

areapunctapercell.soft = [];
areapunctapercell.int = [];
areapunctapercell.stiff = [];

unfilt_sum.soft = [];
unfilt_sum.int = [];
unfilt_sum.stiff = [];

EGFR_sum.soft = [];
EGFR_sum.int = [];
EGFR_sum.stiff = [];

EGFR_percellarea.soft = [];
EGFR_percellarea.int = [];
EGFR_percellarea.stiff = [];

%%
for j = 1:nf
    file = open(files(j).name);
    if contains(files(j).name,'130') == 1
        areapuncta.soft = [areapuncta.soft; [file.data.areapuncta]'];
        areapunctapercell.soft = [areapunctapercell.soft; [file.data.areapuncta_per_areacell]'];
        unfilt_sum.soft = [unfilt_sum.soft; [file.data.unfilt_sum]'];
        EGFR_sum.soft = [EGFR_sum.soft; [file.data.EGFR_sum]'];
        EGFR_percellarea.soft = [EGFR_percellarea.soft; [file.data.EGFR_percell]'];
    elseif contains(files(j).name,'910') == 1
        areapuncta.int = [areapuncta.int; [file.data.areapuncta]'];
        areapunctapercell.int = [areapunctapercell.int; [file.data.areapuncta_per_areacell]'];
        unfilt_sum.int = [unfilt_sum.int; [file.data.unfilt_sum]'];
        EGFR_sum.int = [EGFR_sum.int; [file.data.EGFR_sum]'];
        EGFR_percellarea.int = [EGFR_percellarea.int; [file.data.EGFR_percell]'];
    elseif contains(files(j).name,'4020') == 1
        areapuncta.stiff = [areapuncta.stiff; [file.data.areapuncta]'];
        areapunctapercell.stiff = [areapunctapercell.stiff; [file.data.areapuncta_per_areacell]'];
        unfilt_sum.stiff = [unfilt_sum.stiff; [file.data.unfilt_sum]'];
        EGFR_sum.stiff = [EGFR_sum.stiff; [file.data.EGFR_sum]'];
        EGFR_percellarea.stiff = [EGFR_percellarea.stiff; [file.data.EGFR_percell]'];
    else
        break
    end
end

EGF_per_EGFR.soft = areapuncta.soft ./ EGFR_sum.soft;
EGF_per_EGFR.int = areapuncta.int ./ EGFR_sum.int;
EGF_per_EGFR.stiff = areapuncta.stiff ./ EGFR_sum.stiff;

unfilt_sum_per_EGFR.soft = unfilt_sum.soft ./ EGFR_sum.soft;
unfilt_sum_per_EGFR.int = unfilt_sum.int ./ EGFR_sum.int;
unfilt_sum_per_EGFR.stiff = unfilt_sum.stiff ./ EGFR_sum.stiff;

figure(1)
nbins = 20;
hold on
histogram(areapuncta.soft,nbins)
histogram(areapuncta.int,nbins)
histogram(areapuncta.stiff,nbins)
legend('soft','intermediate','stiff')

  figure(2)
  hold on
  plot(EGFR_sum.stiff,areapuncta.stiff,'o')
  plot(EGFR_sum.int,areapuncta.int,'o')
  plot(EGFR_sum.soft,areapuncta.soft,'o')
  legend('stiff','intermediate','soft')
  
filename = '2021_02_26.mat';
save(filename)