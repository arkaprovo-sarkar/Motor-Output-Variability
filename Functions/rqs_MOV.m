%% Function to perform RQS (paired and unpaired)
% 
% Input:
% data = Target curve across trials (sample points x trials)
function [RQS] = rqs_MOV(data)

t_stat  = zeros(1, size(data, 2));
for trl = 1:size(data, 2)
    mdl = fitlm((data(:,i), mean(data,2))); % inputs (X,Y), column vectors where columns of Y are the response variable while X are the predictior variable
    t_stat(i) = abs(mdl.Coefficients.tStat(2)); % first row is intercept while second row predictor variable
end
            RQS(subject, cond) = mean(t_stat);