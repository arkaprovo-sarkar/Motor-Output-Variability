% Check baseline STD

function [bad_trials] = baseline_STD(data, baseline,times)

%for trl = 1:size(data,2)
bl_idx = find(times>=baseline(1) & times <=baseline(2));
mean_bl = mean(data(bl_idx,:), 1);
mu = mean(mean_bl);
sigma = std(mean_bl);
bad_trials = find (abs(mean_bl - mu) > 3*sigma);
end
