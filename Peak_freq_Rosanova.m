% Natural Freqency: Rosanova 2009 method

% Run Morlet Wavelet TF analysis

clc
clear

subject_cohort = [1 2 11:17 19:27];
condition_cohort = [1 3 6 4];

% Set directory to functions folder
cd D:\ARKO\DATA\MotorOutputVariability\Scripts_Final\Functions

for cond = 1:length(condition_cohort)
    
    % Setting condition dependent baseline (in ms)
    if condition_cohort(cond) == 1 || condition_cohort(cond) == 3
        baseline = [-1250 -750];
    elseif condition_cohort(cond) == 6
        baseline = [-750 -250];
    elseif condition_cohort(cond) == 4
        baseline = [-2249 -1749];
    end

    for subject = subject_cohort
        
        % Load clean EEG data
        eeglab;
        subject_ID = [sprintf('%03d', subject)];
        directory = ['D:\ARKO\DATA\MotorOutputVariability\Preprocessing_Final\cond_',num2str(condition_cohort(cond)),'\'];
        file = [subject_ID, '_EEG_cond', num2str(condition_cohort(cond)),'_v2.set'];
        EEG = pop_loadset('filename',file,'filepath',directory);
        
        % Use custom function for extracting data
        [data_C3, data_M1, LMFP, data_C1, data_C5, data_FC3, data_CP3] = EEG_chselector(EEG);
        
        % Run analysis separately for associated electrodes (check
        % newtimef() for more info on inputs

        % C3
        figure
        [ersp_C3, itc_C3, powbase_C3, times, freqs, erspboot_C3, itcboot, tfdata_C3] = ...
                          newtimef(data_C3, 4000, [-2500 1499], 1000, [3.5], 'baseline', baseline, 'plotphase', 'off', 'padratio', 2, 'freqs', [8 45], 'alpha', 0.05, 'mcorrect', 'fdr', 'naccu', 500,  'timesout', 400);
        % C1
        figure
        [ersp_C1, itc_C1, powbase_C1, times, freqs, erspboot_C1, itcboot, tfdata_C1] = ...
                          newtimef(data_C1, 4000, [-2500 1499], 1000, [3.5], 'baseline', baseline, 'plotphase', 'off', 'padratio', 2, 'freqs', [8 45], 'alpha', 0.05, 'mcorrect', 'fdr', 'naccu', 500,  'timesout', 400);

        % C5
        figure
        [ersp_C5, itc_C5, powbase_C5, times, freqs, erspboot_C5, itcboot, tfdata_C5] = ...
                          newtimef(data_C5, 4000, [-2500 1499], 1000, [3.5], 'baseline', baseline, 'plotphase', 'off', 'padratio', 2, 'freqs', [8 45], 'alpha', 0.05, 'mcorrect', 'fdr', 'naccu', 500,  'timesout', 400);

        % FC3
        figure
        [ersp_FC3, itc_FC3, powbase_FC3, times, freqs, erspboot_FC3, itcboot, tfdata_FC3] = ...
                          newtimef(data_FC3, 4000, [-2500 1499], 1000, [3.5], 'baseline', baseline, 'plotphase', 'off', 'padratio', 2, 'freqs', [8 45], 'alpha', 0.05, 'mcorrect', 'fdr', 'naccu', 500,  'timesout', 400);

        % CP3
        figure
        [ersp_CP3, itc_CP3, powbase_CP3, times, freqs, erspboot_CP3, itcboot, tfdata_CP3] = ...
                          newtimef(data_CP3, 4000, [-2500 1499], 1000, [3.5], 'baseline', baseline, 'plotphase', 'off', 'padratio', 2, 'freqs', [8 45], 'alpha', 0.05, 'mcorrect', 'fdr', 'naccu', 500,  'timesout', 400);

        close all;
        signal_M1(:,:,1) = ersp_C3; signal_M1(:,:,2) = ersp_C1; signal_M1(:,:,3) = ersp_C5; signal_M1(:,:,4) = ersp_FC3; signal_M1(:,:,5) = ersp_CP3;
        ersp_M1 = mean(signal_M1, 3);
        save (['D:\ARKO\DATA\MotorOutputVariability\TF_data\Milano_parameters\cond',num2str(condition_cohort(cond)),'\',num2str(subject_ID),'_cond', num2str(condition_cohort(cond)), '_TF.mat'], 'ersp_M1', 'ersp_C3', 'itc_C3', 'powbase_C3', 'erspboot_C3', 'tfdata_C3', ...
           'ersp_C1', 'itc_C1', 'powbase_C1', 'erspboot_C1', 'tfdata_C1', 'ersp_C5', 'itc_C5', 'powbase_C5', 'erspboot_C5', 'tfdata_C5', 'ersp_FC3', 'itc_FC3', 'powbase_FC3', 'erspboot_FC3', 'tfdata_FC3','ersp_CP3', 'itc_CP3', 'powbase_CP3', 'erspboot_CP3', 'tfdata_CP3', 'times', 'freqs' );
        
       end
end

%%
clc
clear

% Initialize variables
subject_cohort = [1 2 11:17 19:27];
condition_cohort = [1 3 6 4];
nfreq = 38; % No of freq bins from Morlet Wavelet TF analysis
avg_wave= zeros(length(subject_cohort), nfreq, length(condition_cohort));
time_range = [20 200];
%fb_matrix = [8 45; 8 13; 13 30; 30 45];
plot_count = 1; % suplot counter
colour_matrix = [2 197 247; % Rest
                 254 191 4;  % FTV 700
                 252 108 133; % FTV 200
                 18 168 108]; % Movement Onset
bar_colours = colour_matrix./255; % Convert RGB values to range [0, 1]
peak_freq_all = zeros(length(subject_cohort), length(condition_cohort), 2); % Store peak freq and probability of occurence
cond_names = {'Rest', 'FTV700', 'FTV200', "MovtOnset"}
C3_storage = zeros(length(subject_cohort), nfreq,  length(condition_cohort)); % Storing C3 ERSP values (subjects x freqs x conditions)
peak_freq_C3 = zeros(length(subject_cohort),  length(condition_cohort)); % Storing ERSP C3 peaks

% Loop through subjects and conditions
for cond = 1:length(condition_cohort)
figure;

plot_idx = 1;
for subject = 1:length(subject_cohort)

    subject_ID = [sprintf('%03d', subject_cohort(subject))];
    data = open(['D:\ARKO\DATA\MotorOutputVariability\TF_data\Milano_parameters\cond',num2str(condition_cohort(cond)),'\', subject_ID, '_cond', num2str(condition_cohort(cond)), '_TF.mat']);       

    % Get the significant values from bootstrapped ERSP of differnt
    % channels
    erspC3_mask = data.ersp_C3 <= data.erspboot_C3(:,1) | data.ersp_C3 >= data.erspboot_C3(:,2);
    erspC1_mask = data.ersp_C1 <= data.erspboot_C1(:,1) | data.ersp_C1 >= data.erspboot_C1(:,2);
    erspC5_mask = data.ersp_C5 <= data.erspboot_C5(:,1) | data.ersp_C5 >= data.erspboot_C5(:,2);
    erspCP3_mask = data.ersp_CP3 <= data.erspboot_CP3(:,1) | data.ersp_CP3 >= data.erspboot_CP3(:,2);
    erspFC3_mask = data.ersp_FC3 <= data.erspboot_FC3(:,1) | data.ersp_FC3 >= data.erspboot_FC3(:,2);
    
    erspC3_sig = data.ersp_C3.*erspC3_mask;
    erspC1_sig = data.ersp_C1.*erspC1_mask;
    erspC5_sig = data.ersp_C5.*erspC5_mask;
    erspFC3_sig = data.ersp_FC3.*erspFC3_mask;
    erspCP3_sig = data.ersp_CP3.*erspCP3_mask;

    time_window = find(data.times>=time_range(1) & data.times<= time_range(2));

    % Cumulative ERSP over time window
    cpow_C3 = sum(erspC3_sig(:,time_window), 2);
    cpow_C1 = sum(erspC1_sig(:,time_window), 2);
    cpow_C5 = sum(erspC5_sig(:,time_window), 2);
    cpow_FC3 = sum(erspFC3_sig(:,time_window), 2);
    cpow_CP3 = sum(erspCP3_sig(:,time_window), 2);

    C3_storage(subject, :, cond) = cpow_C3; % Store
    peak_freq = data.freqs(find(cpow_C3 == max(cpow_C3)));
    % Handle flat peaks
    if size(peak_freq,2) > 1
        peak_freq = peak_freq(1); % At plateau take the lowest freq value
    end
    % Plotting
    subplot(3, 6, plot_idx)
    plot(data.freqs, cpow_C3, 'LineWidth', 1.5, 'Color',bar_colours(cond,:));
    ylim([min(cpow_C3)-10 max(cpow_C3)+10])
    xlim([4 45])
    xticks([8 15 25 35 45])
    plot_idx = plot_idx + 1;
    hold on;
    plot(peak_freq, max(cpow_C3), 'o', 'MarkerSize', 7);
    ylim([min(cpow_C3)-10 max(cpow_C3)+10])
    peak_freq_C3(subject, cond) = peak_freq;
    title([subject_ID, ' | Cond: ', char(cond_names{cond}), ' | PF: ', num2str(peak_freq)], 'FontSize',12);
end
end

%%
% Visualization
% Bar plot
y_label = 'Hz';
title_name = 'Natural Freq';
var_labels ={'Rest', 'FTV 700', 'FTV 200', 'MovtOnset'};
y_range = [8 45];
p_value = barplot_MOV(peak_freq_C3, title_name, var_labels, y_label, y_range);

