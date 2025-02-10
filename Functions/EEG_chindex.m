%% Channel indexing
%
% Input:
% chanlocs = Data structure loaded from EEGLab EEG.chanlocs
% ch_names = Char array of chanel names needing indexing
%
% Output:
% ch_names = Channel names 
% ch_index = Index of the the requested channel names


function [ch_index, ch_names] = EEG_chidx(chanlocs, ch_names)

ch_index = zeros(size(ch_names));

for i=1:length(ch_names)
    for j = 1:length(chanlocs)
        if strcmp(chanlocs(j).labels , ch_names{i}) == 1
            ch_index(i) = j;
        end
    end
end