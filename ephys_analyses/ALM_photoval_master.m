%% Load appropriate sp_struct

% Load deep or superficial cortical inactivation recordings
%load('ALM_sup_val.mat')
load('ALM_deep_val.mat')

%% Filter sp_struct of neurons that do not meet minimum criteria

params.min_trial_num = 10;
params.min_fr_thresh = 1;

sp_struct_filt = filter_sp_struct_photoval(sp_struct,params);

%% Calculate mean firing rates and smooth

sp_times_filt = {sp_struct_filt.sp_times_filt};
laser_trial = {sp_struct_filt.laser_trial};

% note that, we want to have spare bins to the left and right of what we
% care about in bin_window for smoothing of PSTHs below.
smooth_bin_window = [-0.35 0.65];
real_window = [-0.3 0.6];
fr_window = [0 0.3];
bin_size = 0.01;
filt_window = (smooth_bin_window - real_window)/bin_size;
filt_ind = round([abs(filt_window(1)) + 1; sum(abs(smooth_bin_window))/bin_size - abs(filt_window(1))]);
mean_window = (smooth_bin_window - fr_window)/bin_size;
mean_ind = round([abs(mean_window(1)) + 1; sum(abs(smooth_bin_window))/bin_size - abs(mean_window(1))]);
convert_to_sec = 1;
sp_times_bin = cell(numel(sp_struct_filt), 1);
for i = 1:numel(sp_struct_filt)
    % align and bin firing rates
    sp_times_bin_temp = get_bin_firing_rates(sp_times_filt{i}, [], [], [], smooth_bin_window, bin_size, convert_to_sec);
    sp_times_bin{i}  = sp_times_bin_temp;
end

% create PSTHs by averaging trials over conditions and smooth them
bin_size_back = 3;
bin_size_forward = 0;
mean_firing_rate = nan(numel(sp_struct_filt), 2, size(sp_times_bin{1}, 2) - sum(abs(filt_window)));
mean_firing_rate_win = nan(numel(sp_struct_filt), 2);
for i = 1:numel(sp_times_bin)
    laser_trial_temp = [laser_trial{i}{:}];

    % not the most efficient way to do this...fix later. BSI.
    % min-max
    mean_firing_rate_both = [mean(sp_times_bin{i}(laser_trial_temp == 1, :)) mean(sp_times_bin{i}(laser_trial_temp == 2, :))];
    min_fr = min(mean_firing_rate_both);
    max_fr = max(mean_firing_rate_both);

    for j = 1:2

        % normalize firing rates
        mean_firing_rate_temp = mean(sp_times_bin{i}(laser_trial_temp == (j - 1), :));
        norm_mean_firng_rate_temp = (mean_firing_rate_temp - min_fr)/(max_fr - min_fr);

        % smooth firing rates
        mean_firing_rate_temp = smoothdata(norm_mean_firng_rate_temp, 'movmean', [bin_size_back bin_size_forward]);
        
        % assign firing rates
        mean_firing_rate(i, j, :) = norm_mean_firng_rate_temp(filt_ind(1):filt_ind(2));
        mean_firing_rate_win(i, j) = mean(mean_firing_rate_temp(mean_ind(1):mean_ind(2)));
    end
end

% find neurons that were excited or inhibited
inh = mean_firing_rate_win(:, 2) < mean_firing_rate_win(:, 1);
exc = mean_firing_rate_win(:, 2) > mean_firing_rate_win(:, 1);

% sort mean_firing_rate from highest to lowest mean FR
[~, ind] = sort(mean_firing_rate_win(:, 2), 'ascend');

%% Create heatmaps

figure();
imagesc([squeeze(mean_firing_rate(ind, 2, :))]);
colormap(parula)
clim([0 1])
colorbar

%% Do a hierarchical bootstrap for the IQR and plot for Cortex

sp_times_filt_inh = sp_times_filt(inh);
laser_trial_inh = laser_trial(inh);

num_shuffles = 1000;
mean_firing_rate_win_rand = nan(numel(sp_times_filt_inh), 2, num_shuffles);
num_neurons = numel(sp_times_filt_inh);
parfor x = 1:num_shuffles

    disp(x)

    rand_ind = randi(numel(sp_times_filt_inh), 1, numel(sp_times_filt_inh));
    
    sp_times_filt_rand = sp_times_filt_inh(rand_ind);
    laser_trial_rand = laser_trial_inh(rand_ind);
    
    % note that, we want to have spare bins to the left and right of what we
    % care about in bin_window for smoothing of PSTHs below.
    smooth_bin_window = [-0.35 0.65];
    real_window = [-0.3 0.6];
    fr_window = [0 0.3];
    bin_size = 0.01;
    filt_window = (smooth_bin_window - real_window)/bin_size;
    filt_ind = round([abs(filt_window(1)) + 1; sum(abs(smooth_bin_window))/bin_size - abs(filt_window(1))]);
    mean_window = (smooth_bin_window - fr_window)/bin_size;
    mean_ind = round([abs(mean_window(1)) + 1; sum(abs(smooth_bin_window))/bin_size - abs(mean_window(1))]);
    convert_to_sec = 1;
    sp_times_bin = cell(numel(sp_times_filt_inh), 1);
    for y = 1:num_neurons
        % align and bin firing rates
        sp_times_bin_temp = get_bin_firing_rates(sp_times_filt_rand{y}, [], [], [], smooth_bin_window, bin_size, convert_to_sec);
        sp_times_bin{y}  = sp_times_bin_temp;
    end
    
    % create PSTHs by averaging trials over conditions and smooth them
    bin_size_back = 3;
    bin_size_forward = 0;
    mean_firing_rate = nan(numel(sp_times_filt_inh), 2, size(sp_times_bin{1}, 2) - sum(abs(filt_window)));
    for i = 1:num_neurons
        laser_trial_temp = [laser_trial_rand{i}{:}];

        % not the most efficient way to do this...fix later. BSI.
        % min-max
        mean_firing_rate_both = [mean(sp_times_bin{i}(laser_trial_temp == 1, :)) mean(sp_times_bin{i}(laser_trial_temp == 2, :))];
        min_fr = min(mean_firing_rate_both);
        max_fr = max(mean_firing_rate_both);
        for j = 1:2

            mean_firing_rate_temp = mean(sp_times_bin{i}(laser_trial_temp == (j - 1), :));
            norm_mean_firng_rate_temp = (mean_firing_rate_temp - min_fr)/(max_fr - min_fr);

            % smooth firing rates
            mean_firing_rate_temp = smoothdata(norm_mean_firng_rate_temp, 'movmean', [bin_size_back bin_size_forward]);
            mean_firing_rate_win_rand(i, j, x) = mean(mean_firing_rate_temp(mean_ind(1):mean_ind(2)));
        end
    end
end

med_firing_rate_rand = squeeze(median(mean_firing_rate_win_rand));
iqr_rand = prctile(med_firing_rate_rand', [25 75]);

med_no_laser = median(mean_firing_rate_win(:, 1));
med_laser = median(mean_firing_rate_win(:, 2));

figure();
scatter(1, med_no_laser); hold on;
scatter(2, med_laser); 
errorbar(1, med_no_laser, med_no_laser - iqr_rand(1, 1), iqr_rand(2, 1) - med_no_laser);
errorbar(2, med_laser, med_laser - iqr_rand(1, 2), iqr_rand(2, 2) - med_laser);
xlim([0 3]);
ylim([0 0.5]);