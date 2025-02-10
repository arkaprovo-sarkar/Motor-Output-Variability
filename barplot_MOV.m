% Plots bar graphs and does stats for you

function [p_matrix] = barplot_MOV(plot_data, title_name, var_labels, y_label, y_range, x_label)


    colour_matrix = [2 197 247; % Rest
                     254 191 4;  % FTV 700
                     252 108 133; % FTV 200
                     18 168 108]; % Movement Onset
    bar_colours = colour_matrix./255; % Convert RGB values to range [0, 1]

if length(size(plot_data)) == 2
    %plot_data = plot_data - plot_data(:,1); % For unpaired RQS
    mean_values = mean(plot_data);
    std_errors = std(plot_data)/sqrt(length(plot_data));
    
    % % Box Plot
    % figure;
    % b = boxplot(plot_data, 'Labels', var_labels)
    % ylabel(y_label);
    % title('Boxplot', 'FontSize', 25);
    % xticklabels(var_labels);
    % set(gca, 'FontSize', 10); % Adjust font size of tick labels
    % 
    figure;
    b = bar([1:length(var_labels)], mean_values, 'EdgeColor', 'k', 'LineWidth', 0.7);
    b.FaceColor = "flat";
    b.CData = bar_colours;
    hold on;
    errorbar([1:length(var_labels)], mean_values, std_errors, 'k.', 'LineWidth', 1.5);
    hold off;
    ylabel(y_label,  'FontSize', 18);
    title(title_name, 'FontSize', 18);
    xticklabels(var_labels);
    set(gca, 'FontSize', 18); % Adjust font size of tick labels
    if ~isempty(y_range)
    ylim(y_range);
    end
    
    h_matrix = zeros(size(plot_data,2));
    p_matrix = zeros(size(plot_data,2));


    % Stat matrix (post-hoc)
    for i = 1:size(plot_data,2)
        for j = 1: size(plot_data,2) 
            [p, h] = signrank(plot_data(:,i), plot_data(:,j));
            h_matrix(i, j) = h;
            p_matrix(i, j) = p;
        end
    end

    %Make cloud plot
    for i = 1:size(plot_data, 1)
        hold on;
        scatter(1:size(plot_data,2), plot_data(i,:), 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
        %Connect the points with lines
        %plot(1:size(plot_data,2), plot_data(i,:), 'Color', [0.5 0.5 0.5]);
    end
    % 
elseif length(size(plot_data)) == 3
    
    f_bandnames = {'All (4-45 Hz)', 'Theta (4-8 Hz)', 'Alpha (8-13 Hz)', 'Beta (13-30 Hz)', 'Low Gamma(30-45 Hz)'};
    mean_values = squeeze(mean(plot_data, 1)); % Freq bands x conditions
    std_errors = squeeze(std(plot_data, 1))./sqrt(size(plot_data,1));

    % Plot the bar graph with specified colors
    b = bar(mean_values); % mean values represent avg values x conditions in their dimension
    hold on;
    
    % Specify colours based on conditions
    for i = 1:size(mean_values, 2)
        b(i).FaceColor = bar_colours(i,:);  
    end
    
    % Calculate position the centre of grouped bars and use it to plot error bars
    ngroups = size(mean_values, 1);
    nbars = size(mean_values, 2);
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    for i = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x, mean_values(:,i), std_errors(:,i), '.', 'LineWidth', 1);
    end
    
    % Add legend with correct colors
    labels = {'Rest', 'FTV 700', 'FTV 200', 'MovtOnset'};
    legend(labels, 'Location', 'best', 'FontSize', 18);
    
    % Add labels and title
    xlabel(x_label, 'FontSize', 14);
    ylabel(y_label, 'FontSize', 14);
    title(title_name, 'FontSize', 16);
    set(gca, 'XTick', 1:size(mean_values, 1), 'XTickLabel', f_bandnames, 'FontSize', 18);
    hold off;
    ylim(y_range) % Change based on graphs

end
end