%% EEG time-freq mean
% This function will calculate the mean power across trials within the freq
% window of interest (specified). It will also calculate the average band
% power. It will also do the same with ITC. (Needs single subject data)

% Input:
% data = The time-frequency data information stored by Freq_03_csd.m
% time_range = Vector with first and last element corresponding to the fist
% and last time point respectively (inclusive)
% fb_matrix = Matrix with each row corresponding to the desired freq band with the first and
% second columns corresponding to the first and last freq within the band
% (inclusive)
%
% Output:
% 4 curves (ersp and itc) and 4 band powers/itc avg

function [ersp_curve_C3, ersp_curve_M1, bandpower_C3, bandpower_M1, freqs] = EEG_tfmean(data, time_range, fb_matrix)
    % [ersp_curve_C3, ersp_curve_M1, itc_curve_C3, itc_curve_M1, bandpower_C3, bandpower_M1, avgITC_C3, avgITC_M1, freqs]
    % Get pow-freq curve within specified time-range
    % Target Time
    times_idx = find(data.times >= time_range(1) & data.times <= time_range(2));
    ersp_curve_C3 = mean(data.ersp_C3(:, times_idx), 2);
    ersp_curve_M1 = mean(data.ersp_M1(:, times_idx), 2);
    %itc_curve_C3 = mean(abs(data.itc_C3(:, times_idx)), 2);
    %itc_curve_M1 = mean(abs(data.itc_M1(:, times_idx)), 2);
    freqs = data.freqs;

    % Get average power in different frequency bands
    no_of_bands = size(fb_matrix, 1);
    bandpower_C3 = zeros(1, no_of_bands);
    bandpower_M1 = zeros(1, no_of_bands);
    %avgITC_C3 = zeros(1, no_of_bands);
    %avgITC_M1 = zeros(1, no_of_bands);
    for i=1:no_of_bands
        freq_idx = find(freqs >= fb_matrix(i, 1) & freqs <= fb_matrix(i, 2));
        bandpower_C3(i) = mean(ersp_curve_C3(freq_idx));
        bandpower_M1(i) = mean(ersp_curve_M1(freq_idx));
        %avgITC_C3(i) = mean(itc_curve_C3(freq_idx));
        %avgITC_M1(i) = mean(itc_curve_M1(freq_idx));
    end
end





