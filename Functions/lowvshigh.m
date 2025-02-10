% Low vs High comparisons

% MEP comparisons

clc
clear


subject_cohort = [1 2 11:17 19:27];
condition_cohort = [1 3 6 4];
MEPs_all = zeros(length(subject_cohort), length(condition_cohort), 3); % Last dim for the 3 categories

for cond= 1:length(condition_cohort)

    for subject = 1:length(subject_cohort)

       
        subject_ID = [sprintf('%03d', subject_cohort(subject))];
        
        % Open and sort MEPs
        data = open(['D:\ARKO\DATA\MotorOutputVariability\Preprocessing_Final\cond_', num2str(condition_cohort(cond)), '\', subject_ID, '_MEP_amp_cond', num2str(condition_cohort(cond)), '.mat']);
        MEPs = data.peak_to_peak;
        [sortMEP, sortidx] = sort(MEPs, 'ascend');

        % Select 1/3rd trails (lowerst, mid and highest)
        low_idx = sortidx(1:round(1/3*length(sortidx)));
        high_idx = sortidx(length(sortidx) - (round(1/3*length(sortidx)))+ 1:end);
        mid_idx = setdiff(sortidx, [low_idx; high_idx]); % Removing the high and low trials 

        MEP_low = MEPs(low_idx);
        MEP_mid = MEPs(mid_idx);
        MEP_high = MEPs(high_idx);

        MEPs_all(subject, cond, 1) = mean(MEP_low);
        MEPs_all(subject, cond, 2) = mean(MEP_mid);
        MEPs_all(subject, cond, 3) = mean(MEP_high);

    end
end

MEPs_all = MEPs_all./1000; % conver to mV

% Plot
y_label = 'mV';
title_name = 'MEPs High';
var_labels ={'Rest', 'FTV 700', 'FTV 200', 'MovtOnset'};
y_range = [];
p_value = barplot_MOV(MEPs_all(:,:,3), title_name, var_labels, y_label, y_range);

%% Alpha power comparisons: low med and high

clc
clear

subject_cohort = [1 2 11:17 19:27];
condition_cohort = [1 3 6 4];

win_preTMS = [-105 -5];
win_postTMS = [15 100];
time_matrix = [15 25; 25 35; 40 50; 55 65; 80 120];

fb_matrix = [4 8; 8 13; 13 30; 30 45]
MEPs_all = zeros(length(subject_cohort), length(condition_cohort), 3); % Last dim for the low, med and high groups
zMEPs_all = zeros(length(subject_cohort), length(condition_cohort), 3);
prepow_all = zeros(length(subject_cohort), length(condition_cohort), 3, size(fb_matrix,1));
postpow_all = zeros(length(subject_cohort), length(condition_cohort), 3, size(fb_matrix,1));
TEPs_all = zeros(length(subject_cohort), length(condition_cohort), 6, 3); % 6 peaks, 3 conditions(low mid high)


for cond= 1:length(condition_cohort)

    for subject = 1:length(subject_cohort)

       subject_ID = [sprintf('%03d', subject_cohort(subject))];
       for groups = 1:3
            if groups == 1
                gname = 'low';
            elseif groups == 2
                gname = 'mid';
            elseif groups == 3
                gname = 'high';
            end
           

            % Load Power
            data = open(['D:\ARKO\DATA\MotorOutputVariability\TF_data\Tubingen_parameters\Sorted_trials\cond',num2str(condition_cohort(cond)),'\',num2str(subject_ID),'_cond', num2str(condition_cohort(cond)), '_', gname,'TF.mat']);
            gidx = data.gidx;
            pow_spec = data.ersp_M1;
            freqs = data.freqs;
            times = data.times;
            pretime_idx = find(times>=win_preTMS(1) & times <= win_preTMS(2));
            posttime_idx = find(times>=win_postTMS(1) & times <= win_postTMS(2));

            preersp_avg = mean(pow_spec(:, pretime_idx), 2);
            postersp_avg = mean(pow_spec(:, posttime_idx), 2);

           for f = 1:size(fb_matrix,1)
               f_idx = find(freqs>= fb_matrix(f, 1) & freqs <= fb_matrix(f, 2));
               prepow_all(subject, cond, groups, f) = mean(preersp_avg(f_idx));
               postpow_all(subject, cond, groups, f) = mean(postersp_avg(f_idx));
           end

           % Load MEPs
            data = open(['D:\ARKO\DATA\MotorOutputVariability\Preprocessing_Final\cond_', num2str(condition_cohort(cond)), '\', subject_ID, '_MEP_amp_cond', num2str(condition_cohort(cond)), '.mat']);
            MEPs = data.peak_to_peak;
            MEPs_all(subject, cond, groups) = mean(MEPs(gidx));

            data = open(['D:\ARKO\DATA\MotorOutputVariability\Preprocessing_Final\cond_', num2str(condition_cohort(cond)), '\', subject_ID, '_MEP_zscore_cond', num2str(condition_cohort(cond)), '.mat']);
            zMEPs = data.cond_zscore_MEP;
            zMEPs_all(subject, cond, groups) = mean(zMEPs(gidx));

            % Load TEPs
            data = open(['D:\ARKO\DATA\MotorOutputVariability\TEP_data\cond', num2str(condition_cohort(cond)),'\', subject_ID,'_TEP.mat']);
            M1_signal = data.M1_data;
            t = data.times;
            y = mean(M1_signal(:,gidx), 2);
            TEP_peaks = TEPpeaks_MOV(t', y, time_matrix);
            TEPs_all(subject, cond, :, groups) = TEP_peaks;

       end
    end
end

MEPs_all = MEPs_all./1000;
%% Correlations

figure;
plot_idx = 1;
condition_names = {'Rest', 'FTV 700', 'FTV 200', 'MovtOnset'};
groups = 2; % 1 = low, 2 = mid, 3 = high
f = 2; % 1 = theta, 2= alpha, 3 = beta, 4 = gamma


for cond = 1:length(condition_cohort)
    subplot(1, 4, plot_idx)
    x = prepow_all(:, cond, groups, f);
    y = TEPs_all(:, cond, 6 ,groups);

    %[rho, pval] = corr(x, y, 'Type', 'Spearman');
    [rho, pval] = corr(x, y, 'Type', 'Spearman');
    scatter(x, y, 'filled');
    hold on;
    title_name = [condition_names{cond}];
    title(title_name);
    xlabel('dB');
    ylabel('zscore');
    lsline;
    text(mean(x), mean(y), sprintf('\\rho = %.2f, p = %.4f', rho, pval), 'FontSize', 12, 'Color', 'red');
    hold off;
    plot_idx = plot_idx +1;
end

        