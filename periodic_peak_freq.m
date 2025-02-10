% Periodic Peak Freq:
%
% This script calculates the peak frequency from the fourier power spectrum
% of induced EEG data.

clc
clear

% Initialize variables
subject_cohort = [1 2 11:17 19:27];
condition_cohort = [1 3 6 4];
powpeaks = zeros(length(subject_cohort), length(condition_cohort));
ipowpeaks = zeros(length(subject_cohort), length(condition_cohort));

% Run loop across subjects and conditions
for cond= 1:length(condition_cohort)
    plt_idx= 1; % counter for subplots
    for subject = 1:length(subject_cohort)

        % Load data from EEGLAB
        eeglab;
        subject_ID = [sprintf('%03d', subject_cohort(subject))];
        directory = ['D:\ARKO\DATA\MotorOutputVariability\Preprocessing_Final\cond_',num2str(condition_cohort(cond)),'\'];
        file = [subject_ID, '_EEG_cond', num2str(condition_cohort(cond)),'_v2.set'];
        EEG = pop_loadset('filename',file,'filepath',directory);
        % Store EEG data matrix (channels x samples x trials)
        EEG_signal = EEG.data;
        % Store time vector
        times = EEG.times;
        close; % closing EEGLAB GUI

        % Perform trial wise Fourier Transform of single trial normal and induced activity
        % Initialize variables for Fourier Transform
        Fs = 1000; % Sampling freq
        time_window = [15 500]; % Time window for FFT (in ms)
        time_idx = find(times >= time_window(1) & times <= time_window(2));
        % Identify channel number
        ch_idx = EEG_chindex(EEG.chanlocs, {'C3'}); % Inputs: List of channels, target channel name 
        L = length(time_idx); % Length of signal (no of samples) | Fs/L is the freq resolution
        %n = 2 ^ nextpow2(L); % Add zero padding to increase freq resolution
        n = L; % No zero padding
        f = Fs/L*(0:n/2); % Freq range (0 till Nyquist in inrements of freq res)
        
        % Initialize variables      
        EEG_induced = zeros(size(EEG_signal));
        EEG_erp = mean(EEG_signal, 3); % Mean EEG (across trials)
        fft_signal = zeros(n/2+1, EEG.trials);  fft_induced = zeros(n/2+1, EEG.trials);
        
        % Perform trial-wise FFT on induced signal
        for trl = 1:size(EEG_signal, 3)
            % Induced Signal: EEG(trial) - meanEEG
            EEG_induced(:,:, trl) = EEG_signal(:,:,trl) - EEG_erp;

            % Run fft on induced and normal EEG
            fft_itarget = fft(EEG_induced(ch_idx,time_idx, trl), n);
            fft_target = fft(EEG_signal(ch_idx,time_idx, trl), n);
        
            % Scale amplitude and calculate single side spectrum
            iP2 = abs(fft_itarget/L);
            P2 = abs(fft_target/L);
            iP1 = iP2(1:n/2+1);
            P1 = P2(1:n/2+1);
            iP1(2:end-1) = 2*iP1(2:end-1); % (Info: 1st and last value of P1 correspond to 0 and nyquist freq which don't have any pow in the negative freq space (the other values have mirrored pow so we double them))
            P1(2:end-1) = 2*P1(2:end-1);

            % Store induced and normal power spectrum
            fft_induced(:,trl) = iP1;
            fft_signal(:,trl) = P1;
        end        
        
        % Plot the induced and fooofed spectrum
        subplot(3, 6, plt_idx)   
        f_idx = find(f>=4 & f<= 45); % Freq index vector
        plot(f(f_idx), mean(fft_induced(f_idx, :),2), 'LineWidth',2, 'Color', 'k');
        
        % Initialize variables for FOOOF
        f_idx = find(f>=4 & f<= 45);
        freqs  = f(f_idx);
        % FOOOF is run on mean power spectrum
        pow_mean = mean(fft_signal(f_idx,:), 2);
        ipow_mean = mean(fft_induced(f_idx,:), 2);
        settings = struct();  % Use defaults FOOOF parameters
        f_range = [4, 45];
        % Run FOOOF
        fooof_results = fooof(freqs, pow_mean, f_range, settings, true);
        ifooof_results = fooof(freqs, ipow_mean, f_range, settings, true);
        % Store FOOOFed spectrum (on power spectrum model)
        fooof_plot = (fooof_results.fooofed_spectrum - fooof_results.ap_fit);
        ifooof_plot = (ifooof_results.fooofed_spectrum - ifooof_results.ap_fit);
        % Store FOOOFed spectrum (on raw power spectrum)
        % fooof_plot = (fooof_results.power_spectrum - fooof_results.ap_fit);
        % ifooof_plot = (ifooof_results.power_spectrum - ifooof_results.ap_fit);

        % Visualization
        %hold on
        %plot(f(f_idx), fooof_plot, 'LineWidth',2, 'Color', 'k');
        %hold on;
        %plot(f(f_idx), ifooof_plot, 'LineWidth',2, 'Color', 'r');
        %plot(f(f_idx), ifooof_results.power_spectrum, 'LineWidth',2, 'Color', 'k');
        %hold on;
        %plot(f(f_idx), ifooof_results.fooofed_spectrum, 'LineWidth',2, 'Color', 'r');

        % Find Peaks from FOOOFed spectrum
        subtracted_spec = ifooof_results.fooofed_spectrum - ifooof_results.ap_fit;
        hold on
        plot(f(f_idx), subtracted_spec, 'LineWidth',2, 'Color', 'r');
        [~, max_idx] = max(subtracted_spec);
        ppFreq = freqs(max_idx); % periodic peak frequncies
        % Visualization of periodic peak frequencies
        hold on
        plot(ppFreq, subtracted_spec(max_idx), 'o', 'Color', 'k', 'MarkerSize',13);
        title([subject_ID, ' | PPF: ', num2str(round(ppFreq, 1)),' Hz'])
        xlabel("f (Hz)")
        ylabel("|P1(f)|")
        xlim([4 45])
        p_freq = freqs(find(fooof_plot == max(fooof_plot)));
        ip_freq = freqs(find(ifooof_plot == max(ifooof_plot)));
        % Handle error when peak not present
        if max(fooof_plot) == 0
            p_freq = 0;
        elseif max(ifooof_plot) == 0
            ip_freq = 0;
        end

        % Store FOOOFed peaks for normal and induced specrum
        powpeaks(subject, cond) = p_freq;
        ipowpeaks(subject, cond) = ip_freq;
       
        plt_idx = plt_idx +1;
    end
    close
end



%% Visualization

% Bar plot
y_label = ['Freq (Hz)'];
title_name = 'Peak Freq | Induced | FFT | C3';
var_labels ={'Rest', 'FTV 700', 'FTV 200', 'MovtOnset'};
y_range = [4 45];
p_value = barplot_MOV(ipowpeaks, title_name, var_labels, y_label, y_range);

% Saving images
%export_fig(['PeakFreqs_cond',num2str(condition_cohort(cond)), '.jpg'], '-jpg', '-r300');
%movefile(['PeakFreqs_cond',num2str(condition_cohort(cond)), '.jpg'], target_folder);

% Box Plot
colour_matrix = [2 197 247; % Rest
                     254 191 4;  % FTV 700
                     252 108 133; % FTV 200
                     18 168 108]; % Movement Onset
bar_colours = colour_matrix./255; % Convert RGB values to range [0, 1]
bar_colours = bar_colours * 0.6;

plot_value = [ipowpeaks(:,1); ipowpeaks(:,2); ipowpeaks(:,3); ipowpeaks(:,4)];
group_inx = [ones(1,18), 2.*ones(1,18) 3.*ones(1,18) 4.*ones(1,18)];

h = daboxplot(plot_value,'groups',group_inx, 'xtlabels', var_labels, 'fill', 0, 'scatter',2, 'scattercolors',{'w', 'k'}, 'colors',bar_colours, 'scattersize', 30, 'mean', 1, ...
    'boxwidth', 2)
xlim([0 5])
set(h.mn,'LineWidth',3); % customize mean lines
set(h.bx,  'LineWidth',2);
set(h.sc, 'LineWidth',1);

ylabel('Hz', 'FontSize', 14);
set(gca,'FontSize',14);
ylim([4 45]);
%title ('Periodic Peak Frequency | C3', 'FontSize',16);

fig.Position = [100, 100, 1200, 600]; % [a, b, c, d]; a = no. of pixels from left, b = no. of pexels from bottom, c = width of fig, d = height of fig]

%save_path = 'D:\ARKO\DATA\MotorOutputVariability\Final_Images\PeakFreq_FFT_FOOOF.jpeg';
%exportgraphics(gcf, save_path, 'BackgroundColor', 'white', 'Resolution', 660);

