% Peak frequency analysis

clc
clear

subject_cohort = [1 2 11:17 19:27];
condition_cohort = [1 3 6 4];
nfreq = 42; % No of freq bins
avg_wave= zeros(length(subject_cohort), nfreq, length(condition_cohort));

time_range = [15 140];
%fb_matrix = [8 45; 8 13; 13 30; 30 45];
plot_count = 1;

colour_matrix = [2 197 247; % Rest
                 254 191 4;  % FTV 700
                 252 108 133; % FTV 200
                 18 168 108]; % Movement Onset
bar_colours = colour_matrix./255; % Convert RGB values to range [0, 1]

peak_freq_all = zeros(length(subject_cohort), length(condition_cohort), 2); % Store peak freq and probability of occurence

freqs_storage_rest = [];
freqs_storage_ftv700 = [];
freqs_storage_ftv200 = [];
freqs_storage_movtonset = [];


for cond = 1:length(condition_cohort)
figure;
if condition_cohort(cond) == 1 || condition_cohort(cond) == 3
    baseline = [-1250 -750];
    epoch_window = [-1500 1000];
elseif condition_cohort(cond) == 6
    baseline = [-750 -250];
    epoch_window = [-1500 1000];
elseif condition_cohort(cond) == 4
    baseline = [-2249 -1749];
    epoch_window = [-2500 1000];
end
plot_count = 1;

for subject = 1:length(subject_cohort)

    subject_ID = [sprintf('%03d', subject_cohort(subject))];
    data = open(['D:\ARKO\DATA\MotorOutputVariability\TF_data\Tubingen_parameters\cond',num2str(condition_cohort(cond)),'\', subject_ID, '_cond', num2str(condition_cohort(cond)), '_TF.mat']);       
    
    % trial adjustment
    window = find(data.times>=epoch_window(1) & data.times <= epoch_window(2));
    data.tfdata_C3 = data.tfdata_C3(:,window,:);
    data.tfdata_C1 = data.tfdata_C1(:,window,:);
    data.tfdata_C5 = data.tfdata_C5(:,window,:);
    data.tfdata_FC3 = data.tfdata_FC3(:,window,:);
    data.tfdata_CP3 = data.tfdata_CP3(:,window,:);
    data.times = data.times(window);

    
    alltfX(:,:,:,1) = data.tfdata_C3; alltfX(:,:,:,2) = data.tfdata_C1; alltfX(:,:,:,3) = data.tfdata_C5; alltfX(:,:,:,4) = data.tfdata_CP3; alltfX(:,:,:,5) = data.tfdata_FC3;
    freqs = data.freqs;
    
    %From 8 Hz onwards (38 BINS, OTHERWISE 42)
    % alltfX = alltfX(5:end,:,:,:);
    % freqs = freqs(5:end);


    stbl_all_elec = zeros(size(alltfX,1), size(alltfX,2), size(alltfX,3), 5); % Initializing storage variable 
        
    for i =1:5 % Check through the 5 electrodes
        
        % First convert complex values to power values P
        P = alltfX(:,:,:,i).*conj(alltfX(:,:,:,i));
        P_norm = zeros(size(P));
        P_norm_bl = zeros(size(P));
        bl_P = zeros(size(P,3),nfreq);
         % Find baseline window
         bl_idx = find(data.times >= baseline(1) & data.times <= baseline(2))
        % Single Trial normalization on power values
        for j = 1:size(P,3) % Looping through each trial
            
            for freq_trl = 1:size(P,1) % Looping through each freq value
                pow_row = P(freq_trl, :, j);
                mean_P = mean(pow_row);
                P_norm(freq_trl, :, j) = (pow_row)./mean_P;
                %P_norm(freq_trl, :, j) = pow_row;

               
                bl_P(j, freq_trl) = mean(P_norm(freq_trl,bl_idx,j));

                % Subtract baseline values within each trial
                P_norm_bl(freq_trl, :, j) = P_norm(freq_trl, :, j)./ bl_P(j, freq_trl);
            end
        end
        % Time for subtracting baseline avg from trials
        for j = 1:size(P,3) % Looping through each trial
            % for freq_trl = 1:size(P,1) % Looping through each freq value
            %     st_P = P_norm(freq_trl,:,j); % Trial wise value
            %     P_norm_bl(freq_trl,:,j) = st_P./ bl_P(j, freq_trl); % Baseline corrected (divided by baseline values)
            % end
            % Now save the curves for the different frequencies for
            % every electrode within a condition (log scale)
            stbl_all_elec(:,:,j, i) = 10*log10(P_norm_bl(:,:,j));
        end
    end
    %stbl_all = mean(stbl_all_elec, 4); % Average the frequency across all the electrodes
    stbl_all = stbl_all_elec(:,:,:,1); % Only C3
    
    target_idx = find(data.times >= time_range(1) & data.times <= time_range(2));
    stbl_ersp_all = squeeze(mean(stbl_all(:, target_idx, :), 2));
    
    peak_freq = zeros(size(stbl_ersp_all, 2),1);
    peak_freq_pow = zeros(size(stbl_ersp_all, 2),1);
    
    cumP(subject, :, cond) = sum(stbl_ersp_all,2); % cumulative waveform
    avg_wave(subject,:, cond) = mean(stbl_ersp_all,2); % Average waveform (to get power at peak freq)

    for trl = 1:size(stbl_ersp_all, 2) % Loop through trials
        peak_freq_pow(trl) = max(stbl_ersp_all(:,trl));
        peak_freq(trl) = freqs(find(stbl_ersp_all(:, trl) == peak_freq_pow(trl)));
    end

    % Plot average waveform and peak freq across trials
    
    subplot(3,6,plot_count)
    %

    h = histfit(peak_freq, 9);
    h(1).FaceColor = bar_colours(cond,:);
    h(2).Color = bar_colours(cond,:)*0.5;
    xticks([8 15 25 35 45]);
    %counts = h.Values;
    %binEdges = h.BinEdges;
    normalpd = fitdist(peak_freq, 'Lognormal');
    hold on
    
    %[~, maxIndex] = max(counts);
    %peak_freq_bin = [binEdges(maxIndex), binEdges(maxIndex + 1)];
    %prob_peak = round(max(counts)/sum(counts)*100,2);
    %title([subject_ID, ' | PF: ', num2str(normalpd.mu),' Hz | Std: ', num2str(normalpd.std)]);
    title([subject_ID, ' | PF: ', num2str(median(peak_freq)),' Hz']);
    %ylim([0 55]);
    xlim([8 45]);
    xlabel('Max Frequency (Hz)');
    ylabel('No. of Trials')
    %xline(mean(peak_freq_bin), 'Color', 'k', 'LineWidth', 2, 'LineStyle','--');
    %xline(normalpd.mu, 'Color', 'k', 'LineWidth', 2, 'LineStyle','--')
    xline(median(peak_freq), 'Color', 'k', 'LineWidth', 2, 'LineStyle','--')
    
    plot_count= plot_count + 1;
    %close all;


    peak_freq_all(subject, cond, 1) = median(peak_freq);
    % peak_freq_all(subject, cond, 2) = prob_peak;
    %peak_freq_all(subject, cond, 1) = normalpd.mu; % Mean 
    peak_freq_all(subject, cond, 2) = normalpd.std; % STD
    if cond == 1
        freqs_storage_rest{subject} = peak_freq;
    elseif cond == 2
        freqs_storage_ftv700{subject} = peak_freq;
    elseif cond == 3
        freqs_storage_ftv200{subject} = peak_freq;
    elseif cond == 4
        freqs_storage_movtonset{subject} = peak_freq;
    end
    clear alltfX
end

save_path = 'D:\ARKO\DATA\MotorOutputVariability\Final_Images\TF_bandpower\Peak_freq_hist_MovtOnset.jpeg';
% exportgraphics(gcf, save_path, 'BackgroundColor', 'white', 'Resolution', 660);
end % cond end

%% Look at avg wave
peak_freq_cum = zeros(18, 4);
natf_power = zeros(18,4);
for cond = 1:4
    figure;
for i = 1: 18
    subplot(3,6,i)
    plot(freqs, cumP(i,:,cond))
    %plot(freqs, avg_wave(i,:,cond))
    %ylim([-1.5 2.5]);
    xlim([4 45]);
    [~, idx] = max(cumP(i,:,cond));
    peak_freq_cum(i, cond) = freqs(idx);
    natf_power(i, cond) = avg_wave(i,idx,cond);
end
end
   
%%


%save('Peak_freqs_allsub.mat', 'peak_freq_all', 'freqs_storage_rest', 'freqs_storage_ftv700', 'freqs_storage_ftv200', 'freqs_storage_movtonset');


colour_matrix = [2 197 247; % Rest
                 254 191 4;  % FTV 700
                 252 108 133; % FTV 200
                 18 168 108]; % Movement Onset
bar_colours = colour_matrix./255; % Convert RGB values to range [0, 1]

% Bar plot
y_label = 'Hz';
title_name = 'Peak Freq';
var_labels ={'Rest', 'FTV 700', 'FTV 200', 'MovtOnset'};
y_range = [];
p_value = barplot_MOV(peak_freq_cum, title_name, var_labels, y_label, y_range);

%% Violin Plot
plot_value = [peak_freq_cum(:,1); peak_freq_cum(:,2);peak_freq_cum(:,3);peak_freq_cum(:,4)];
group_inx = [ones(1,18), 2.*ones(1,18) 3.*ones(1,18) 4.*ones(1,18)];

h = daviolinplot(plot_value,'groups',group_inx, 'xtlabels', var_labels, 'outsymbol','k+','scatter',2,'jitter',1,...
    'box',1,'color',bar_colours, 'boxcolors','same','scattercolors','same',...
    'boxspacing',1.1);

ylabel('Peak Frequency', 'FontSize', 18);
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.2]); % make more space for the legend
set(gca,'FontSize',18);
ylim([0 45]);
ax = gca;
% xtickangle(ax, 45);
fig = gcf;
fig.Position = [100, 100, 1200, 600]; % [a, b, c, d]; a = no. of pixels from left, b = no. of pexels from bottom, c = width of fig, d = height of fig]

% save_path = 'D:\ARKO\DATA\MotorOutputVariability\Final_Images\TF_bandpower\Peak_freq_all.jpeg';
%save_path = ['D:\ARKO\DATA\ICCN_poster\peakFrequency.jpeg'];
%exportgraphics(gcf, save_path, 'BackgroundColor', 'white', 'Resolution', 660);

[p,tbl,stats] = friedman(peak_freq_cum,1)


p = vartestn(squeeze(peak_freq_all(:,:,1)),'TestType','LeveneAbsolute');

% Check the mean
pf_means = median(peak_freq_all(:,:,1),1)

%% DO testing for individual subjects

for i= 1:18
    [h_ftv200(i), p_ftv200(i)] = kstest2(freqs_storage_ftv200{i}, freqs_storage_rest{i});
    [h_ftv700(i), p_ftv700(i)] = kstest2(freqs_storage_ftv700{i}, freqs_storage_rest{i});
    [h_movtonset(i), p_movtonset(i)] = kstest2(freqs_storage_movtonset{i}, freqs_storage_rest{i});
end

%% FOOOF on ERSP

clc
clear

subject_cohort = [1 2 11:17 19:27];
condition_cohort = [1 3 6 4];
nfreq = 42; % No of freq bins
avg_wave= zeros(length(subject_cohort), nfreq, length(condition_cohort));

time_range = [15 100];
%fb_matrix = [8 45; 8 13; 13 30; 30 45];
plot_count = 1;

colour_matrix = [2 197 247; % Rest
                 254 191 4;  % FTV 700
                 252 108 133; % FTV 200
                 18 168 108]; % Movement Onset
bar_colours = colour_matrix./255; % Convert RGB values to range [0, 1]

peak_freq_all = zeros(length(subject_cohort), length(condition_cohort), 2); % Store peak freq and probability of occurence

freqs_storage_rest = [];
freqs_storage_ftv700 = [];
freqs_storage_ftv200 = [];
freqs_storage_movtonset = [];
bl_P_all = [];

ppf_all = zeros(length(subject_cohort), length(condition_cohort));
ppf_pow_all = zeros(size(ppf_all));


%for cond = 1:length(condition_cohort)
cond = 1
figure;
if condition_cohort(cond) == 1 || condition_cohort(cond) == 3
    baseline = [-1250 -750];
    epoch_window = [-1500 1000];
elseif condition_cohort(cond) == 6
    baseline = [-750 -250];
    epoch_window = [-1500 1000];
elseif condition_cohort(cond) == 4
    baseline = [-2249 -1749];
    epoch_window = [-2500 1000];
end
plot_count = 1;


for subject = 1:length(subject_cohort)

    subject_ID = [sprintf('%03d', subject_cohort(subject))];
    data = open(['D:\ARKO\DATA\MotorOutputVariability\TF_data\Tubingen_parameters\cond',num2str(condition_cohort(cond)),'\', subject_ID, '_cond', num2str(condition_cohort(cond)), '_TF.mat']);       
    
    % trial adjustment
    window = find(data.times>=epoch_window(1) & data.times <= epoch_window(2));
    data.tfdata_C3 = data.tfdata_C3(:,window,:);
    data.tfdata_C1 = data.tfdata_C1(:,window,:);
    data.tfdata_C5 = data.tfdata_C5(:,window,:);
    data.tfdata_FC3 = data.tfdata_FC3(:,window,:);
    data.tfdata_CP3 = data.tfdata_CP3(:,window,:);
    data.times = data.times(window);

    
    alltfX(:,:,:,1) = data.tfdata_C3; alltfX(:,:,:,2) = data.tfdata_C1; alltfX(:,:,:,3) = data.tfdata_C5; alltfX(:,:,:,4) = data.tfdata_CP3; alltfX(:,:,:,5) = data.tfdata_FC3;
    freqs = data.freqs;

    stbl_all_elec = zeros(size(alltfX,1), size(alltfX,2), size(alltfX,3), 5); % Initializing storage variable 
        
    for i =1:5 % Check through the 5 electrodes
    %i = 1;
        
        % First convert complex values to power values P
        P = alltfX(:,:,:,i).*conj(alltfX(:,:,:,i));
        P_norm = zeros(size(P));
        P_norm_bl = zeros(size(P));
        bl_P = zeros(size(P,3),nfreq, 5);
         % Find baseline window
         bl_idx = find(data.times >= baseline(1) & data.times <= baseline(2))
        % Single Trial normalization on power values
        for j = 1:size(P,3) % Looping through each trial
            
            for freq_trl = 1:size(P,1) % Looping through each freq value
                pow_row = P(freq_trl, :, j);
                mean_P = mean(pow_row);
                %P_norm(freq_trl, :, j) = (pow_row)./mean_P; % With
                %normalization
                P_norm(freq_trl, :, j) = pow_row; % Without normalization

               
                bl_P(j, freq_trl, i) = mean(P_norm(freq_trl,bl_idx,j));

                % Subtract baseline values within each trial
                P_norm_bl(freq_trl, :, j) = P_norm(freq_trl, :, j)./ bl_P(j, freq_trl, i);
            end
        end
        % Time for subtracting baseline avg from trials
        for j = 1:size(P,3) % Looping through each trial
            % for freq_trl = 1:size(P,1) % Looping through each freq value
            %     st_P = P_norm(freq_trl,:,j); % Trial wise value
            %     P_norm_bl(freq_trl,:,j) = st_P./ bl_P(j, freq_trl); % Baseline corrected (divided by baseline values)
            % end
            % Now save the curves for the different frequencies for
            % every electrode within a condition (log scale)
            % stbl_all_elec(:,:,j, i) = 10*log10(P_norm_bl(:,:,j)); % Saved in log scale
            stbl_all_elec(:,:,j, i) = P_norm_bl(:,:,j); % Saved in power values
        end
    end % Electrode loop

    stbl_all = mean(stbl_all_elec, 4); % Average the frequency across all the electrodes
    %stbl_all = stbl_all_elec(:,:,:,1); % Only C3
    
    target_idx = find(data.times >= time_range(1) & data.times <= time_range(2));
    stbl_ersp_all = squeeze(mean(stbl_all(:, target_idx, :), 2));
    
    freqs = freqs'; % freqs should be in rows (same as PSD)
    % Obtain peak freq (FOOOF)
    psd = mean(stbl_ersp_all,2);
    settings = struct();  % Use defaults
    f_range = [4, 45];
    fooof_results = fooof(freqs, psd, f_range, settings, true);
    AP_comp = fooof_results.ap_fit;
    FOOOF_spec = fooof_results.fooofed_spectrum;
    Pow_spec = fooof_results.power_spectrum;
    fooof_peaks = fooof_results.peak_params(:,:);
    
    subplot(3,6,plot_count)
    standard_error = std(stbl_ersp_all')./size(stbl_ersp_all,2);
    % Shaded error bars
    plot(freqs, Pow_spec, 'Color',bar_colours(cond,:), 'LineWidth',1.5);    
    hold on;
    lineProps = {'-','color',bar_colours(cond,:),'markerfacecolor',bar_colours(cond,:)};
    shadedErrorBar(freqs, Pow_spec, standard_error, lineProps , 1);
    hold on;
    plot(freqs, FOOOF_spec, 'LineWidth', 1.5, 'Color','r'); hold on;
    plot(freqs, AP_comp, 'LineStyle','--', 'LineWidth',1.5, 'Color', 'k');
    plot(freqs, (FOOOF_spec-AP_comp), 'LineWidth',1.5, 'Color', 'r');

    % Get Peak freq
    
    if ~isempty(fooof_peaks)

        freq_peaks = fooof_peaks(:,1);
        ppFreq = freq_peaks(find(fooof_peaks(:,2) == max(fooof_peaks(:,2))));
        max_pow = max(fooof_peaks(:,2));
    else
        max_pow = max(Pow_spec-AP_comp);
        ppFreq = freqs(find((Pow_spec-AP_comp) == max_pow));
    end

    plot(ppFreq, max_pow, 'o', 'Color', 'k', 'MarkerSize',13);
    title([subject_ID, ' | PPF: ', num2str(ppFreq),' Hz']);
    %ylim([0 55]);
    xlim([4 45]);
    xlabel('Frequency (Hz)');
    ylabel('Log Power')    
    plot_count= plot_count + 1;
    %close all;
    % Store stuff
    ppf_all(subject,cond) = ppFreq;
    ppf_pow_all(subject, cond) = max_pow;
 
    clear alltfX
end % Subject end
%end % Condition end

%save_path = 'D:\ARKO\DATA\MotorOutputVariability\Final_Images\FOOOF\MovtOnset_M1.jpeg';
%exportgraphics(gcf, save_path, 'BackgroundColor', 'white', 'Resolution', 660);

%% PPF measure

% Bar plot (Freq)
y_label = 'Hz';
title_name = 'Periodic Peak Freq';
var_labels ={'Rest', 'FTV 700', 'FTV 200', 'MovtOnset'};
y_range = [];
p_value = barplot_MOV(ppf_all, title_name, var_labels, y_label, y_range);

[p,tbl,stats] = friedman(ppf_all,1)

% Bar plot (Power)
y_label = 'Log Power';
title_name = 'Periodic Peak Freq Power';
var_labels ={'Rest', 'FTV 700', 'FTV 200', 'MovtOnset'};
y_range = [];
p_value = barplot_MOV(ppf_pow_all, title_name, var_labels, y_label, y_range);

[p,tbl,stats] = friedman(ppf_pow_all,1)




       
