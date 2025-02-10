% Calculate average band power


function [bandpower_M1_all, avgITC_M1_all, ersp_curve_M1_all, itc_curve_M1_all] = bandpow(subject_cohort, condition_cohort, time_range, fb_matrix)
cond_name = {'Rest', 'FTV 700', 'FTV 200', 'MovtOnset'};
% Initialize storage variables (there are 42 discrete frequencies)
ersp_curve_M1_all = zeros(length(subject_cohort), 42, length(condition_cohort));
itc_curve_M1_all = zeros(length(subject_cohort), 42, length(condition_cohort));

bandpower_M1_all = zeros(length(subject_cohort), size(fb_matrix, 1), length(condition_cohort));
avgITC_M1_all = zeros(length(subject_cohort), size(fb_matrix, 1), length(condition_cohort));

avg_TF_all = zeros(42, 400, 18);

for cond = 1:length(condition_cohort)
    plot_idx = 1;
    for subject = 1:length(subject_cohort)
        subject_ID = [sprintf('%03d', subject_cohort(subject))];
        data = open(['D:\ARKO\DATA\MotorOutputVariability\TF_data\Tubingen_parameters\cond',num2str(condition_cohort(cond)),'\', subject_ID, '_cond', num2str(condition_cohort(cond)), '_TF.mat']);       
        [ersp_curve_C3, ersp_curve_M1, itc_curve_C3, itc_curve_M1, bandpower_C3, bandpower_M1, avgITC_C3, avgITC_M1, freqs] = EEG_tfmean(data, time_range, fb_matrix);      
        % Store values
        ersp_curve_M1_all(subject, :, cond) = ersp_curve_M1';
        itc_curve_M1_all(subject, :, cond) = itc_curve_M1';       
        bandpower_M1_all(subject, :, cond) = bandpower_M1;
        avgITC_M1_all(subject, :, cond) = avgITC_M1;
        
        avg_TF_all(:, :, subject) = data.ersp_M1;
        % Visualize data
        % subplot(4,5,plot_idx)
        % imagesc(data.times, data.freqs, data.ersp_M1)
        % title(['MOV', subject_ID,' | Cond:', num2str(condition_cohort(cond))]);
        % set(gca, 'YDir', 'normal', 'Box', 'off'); % Flip the y-axis and remove box outline
        % colormap('jet');
        % colorbar;
        % x_point = 0; % Time in sec
        % y_limits = ylim; % Get the y-axis limits
        % xlim ([-500 500]);
        % line([x_point x_point], y_limits, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
        % line([(x_point-50) (x_point-50)], y_limits, 'Color', 'k', 'LineWidth', 1.5);
        % % cset colourbar to 0
        % caxis([-max(clim) max(clim)]);
        plot_idx = plot_idx + 1;    
    end
    % % Plot average TF
    % imagesc(data.times, data.freqs, squeeze(mean(avg_TF_all, 3)));
    % title(['Average TF | Cond: ', cond_name(cond)]);
    % set(gca, 'YDir', 'normal', 'Box', 'off'); % Flip the y-axis and remove box outline
    % colormap('jet');
    % colorbar;
    % x_point = 0; % Time in sec
    % y_limits = ylim; % Get the y-axis limits
    % xlim ([-500 500]);
    % line([x_point x_point], y_limits, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
    % line([(x_point-50) (x_point-50)], y_limits, 'Color', 'k', 'LineWidth', 1.5);
    % % cset colourbar to 0
    % caxis([-6.2556 6.2556]);
end

end    
