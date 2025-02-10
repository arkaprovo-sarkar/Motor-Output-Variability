% Band power t-test
% Compute post hoc t test on band power data

function [] = bandpow_testing(band_powers)

% Post hoc-t test
h_matrix = zeros(size(band_powers,3));
p_matrix = zeros(size(band_powers,3));
f_bandnames = {'All (4-45 Hz)', 'Theta (4-8 Hz)', 'Alpha (8-13 Hz)', 'Beta (13-30 Hz)', 'Low Gamma(30-45 Hz)'};    
target_data = band_powers;
for fband = 1:size(target_data, 2)
    for i = 1:size(target_data, 3)
        for j = 1:size(band_powers,3)
            x = squeeze(target_data(:,fband,i));
            y = squeeze(target_data(:,fband,j));
            bf_alpha = 0.05/((size(band_powers,3)*(size(band_powers,3) - 1))/2);
            [p, h] = signrank( x , y, 'alpha', bf_alpha);   
            %[h, p] = ttest( x , y, 'alpha', bf_alpha); 
            h_matrix(i, j) = h;
            if p <= 0.05 % Putting all significant values less than 0.05
                p_matrix(i, j) = round(p, 4);
            else
                p_matrix(i, j) = 1;
            end
        end
    end
    
    % Visualization
    figure;
    imagesc(p_matrix);
    
    custom_colormap = [
        0.9 0.9 0.9; % Light grey color for near significance
        1, 1, 1;  % White colour for non-significant values    
    ];
    colormap(custom_colormap);
    
    for i = 1:length(p_matrix)
        for j = 1:length(p_matrix)
            if p_matrix(i, j) <= bf_alpha
            text(j, i, num2str(p_matrix(i, j)), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'Color', 'b', 'FontSize', 14);
            elseif p_matrix(i, j) > bf_alpha && p_matrix(i, j) <= 0.05
                text(j, i, num2str(p_matrix(i, j)), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'Color', 'k', 'FontSize', 10);
            end
        end
    end
    
    title(['ERSP Band Power ', '| Freq Band: ', char(f_bandnames(fband))], 'FontSize', 12);
    xticklabels({'', 'Rest','', 'FTV 700', '', 'FTV 200', '', 'MovtOnset'});
    yticklabels({'', 'Rest','', 'FTV 700', '', 'FTV 200', '', 'MovtOnset'});
end
end