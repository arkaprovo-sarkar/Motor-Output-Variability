%% Channel selector function
% This function returns the C3 and M1 timeseries data.
%
% Input:
% EEG = Data structure loaded from EEGLab
%
% Output:
% data_C3 = Timeseries data from C3 electrode
% data_M1 = Timeseries data from average of electrodes: C3, C1, C5, FC3,
% CP3

function [data_C3, data_M1, LMFP_data, C1_signal, C5_signal, FC3_signal, CP3_signal] = EEG_chselector(EEG)

 % Identify channels and extract data
        for i=1:length(EEG.chanlocs)

            if strcmp( EEG.chanlocs(i).labels , 'C3') == 1
                channel_C3 = i;
            end
            if strcmp( EEG.chanlocs(i).labels , 'FC1') == 1
                channel_FC1 = i;
            end
            if strcmp( EEG.chanlocs(i).labels , 'CP1') == 1
                channel_CP1 = i;
            end            
            if strcmp( EEG.chanlocs(i).labels , 'FC5') == 1
                channel_FC5 = i;
            end            
            if strcmp( EEG.chanlocs(i).labels , 'CP5') == 1
                channel_CP5 = i;
            end
            if strcmp( EEG.chanlocs(i).labels , 'C1') == 1
                channel_C1 = i;
            end  
            if strcmp( EEG.chanlocs(i).labels , 'C5') == 1
                channel_C5 = i;
            end
            if strcmp( EEG.chanlocs(i).labels , 'CP3') == 1
                channel_CP3 = i;
            end
            if strcmp( EEG.chanlocs(i).labels , 'FC3') == 1
                channel_FC3 = i;
            end              
        end

        C3_signal = EEG.data(channel_C3,:,:);
        FC1_signal = EEG.data(channel_FC1,:,:);
        CP1_signal = EEG.data(channel_CP1,:,:);
        FC5_signal = EEG.data(channel_FC5,:,:);
        CP5_signal = EEG.data(channel_CP5,:,:);
        C1_signal = EEG.data(channel_C1,:,:);
        C5_signal = EEG.data(channel_C5,:,:);
        CP3_signal = EEG.data(channel_CP3,:,:);
        FC3_signal = EEG.data(channel_FC3,:,:);

        % Data to be fed to newtimef()
        data_C3 = C3_signal;
        data_M1 = (C3_signal + C1_signal + C5_signal + CP3_signal + FC3_signal)./5;
        
        LMFP_data = sqrt((C3_signal.^2 + C1_signal.^2 + C5_signal.^2 + CP3_signal.^2 + FC3_signal.^2)./5);
end