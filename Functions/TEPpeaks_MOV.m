% Calculate peaks
% x = times
% y = M1 or C3 mean data
function [TEP_mat] = TEPpeaks_MOV(x, y, time_matrix)

        time = x;
        % Find peaks in the data
        derivative_signal = diff(y) ./ diff(x); % Derivative of signal
        change_points = find(diff(sign(derivative_signal)) ~= 0) + 1; % Points where derivative (slope) changes sign
        time_of_changes = x(change_points);
        peak_index = find(ismember(x, time_of_changes));
        peak_index = peak_index(x(peak_index) > 10 & x(peak_index) < 130);
        
        % Input only the first peak
        N15_flag = false;
        P30_flag = false;
        N45_flag = false;
        P60_flag = false;
        N100_flag = false;
             
        % Calculate the peaks
        % Check this loop.
        for i = 1:length(peak_index)
            peak_saved = false;
            if time(peak_index(i)) > time_matrix(1,1) && time(peak_index(i)) < time_matrix(1,2) && N15_flag == false && peak_saved == false && derivative_signal(peak_index(i)-1) < 0  % makes it a trough
                N15 =  y(peak_index(i));
               
                N15_flag = true;
                peak_saved = true;
            elseif N15_flag == false
                target_time = find(x >= time_matrix(1,1) & x<= time_matrix(1,2));
                N15 =  mean(y(target_time));
            end
            if time(peak_index(i)) > time_matrix(2,1) && time(peak_index(i)) < time_matrix(2,2) && P30_flag == false && peak_saved == false && derivative_signal(peak_index(i)-1) > 0  % makes it a peak
                P30 =  y(peak_index(i));
               
                P30_flag = true;
                peak_saved = true;
            elseif P30_flag == false
                target_time = find(x >= time_matrix(2,1) & x<= time_matrix(2,2));
                P30 =  mean(y(target_time));
               
            end
            if time(peak_index(i)) > time_matrix(3,1) && time(peak_index(i)) < time_matrix(3,2) && N45_flag == false && peak_saved == false && derivative_signal(peak_index(i)-1) < 0  % makes it a trough
                N45 =  y(peak_index(i));
               
                N45_flag = true;
                peak_saved = true;
            elseif N45_flag == false
                target_time = find(x >= time_matrix(3,1) & x<= time_matrix(3,2));
                N45 =  mean(y(target_time));
               
            end
            if time(peak_index(i)) > time_matrix(4,1) && time(peak_index(i)) < time_matrix(4,2) && P60_flag == false && peak_saved == false && derivative_signal(peak_index(i)-1) > 0  % makes it a peak
                P60 =  y(peak_index(i));
                
                P60_flag = true;
                peak_saved = true;
            elseif P60_flag == false
                target_time = find(x >= time_matrix(4,1) & x<= time_matrix(4,2));
                P60 =  mean(y(target_time));
             
            end
            if time(peak_index(i)) > time_matrix(5,1) && time(peak_index(i)) < time_matrix(5,2) && N100_flag == false && peak_saved == false && derivative_signal(peak_index(i)-1) < 0  % makes it a trough
                N100 =  y(peak_index(i));
                
                N100_flag = true;
                peak_saved = true;
            elseif N100_flag == false
                target_time = find(x >= time_matrix(5,1) & x<= time_matrix(5,2));
                N100 =  mean(y(target_time));
               
            end  
        end
    
        % figure;
        % plot(x, y);
        % hold on;
        % plot(time(peak_index), y(peak_index), 'ro', 'MarkerSize', 10);
        % hold off;
        % xlabel('x');
        % ylabel('y');
        % title(['Subject ',num2str(subject)]);

        N15P30 = P30 - N15;

        TEP_mat = [N15; P30; N45; P60; N100; N15P30];
 
end


