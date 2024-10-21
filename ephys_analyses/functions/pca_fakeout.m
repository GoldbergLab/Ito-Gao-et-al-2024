function pca_fakeout(sp_vars_struct, lick_num, num_shuffles)

sp_times_filt = {sp_vars_struct.sp_times_filt};
sp_contact_filt = {sp_vars_struct.sp_contact_filt};
fakeout_trial_filt = {sp_vars_struct.fakeout_trial_filt};
prot_filt = {sp_vars_struct.prot_filt};

% note that, we want to have spare bins to the left and right of what we
% care about in bin_window for smoothing of PSTHs below.
align_lick_num = 2;
smooth_bin_window = [-0.35 0.35];
pca_bin_window = [-0.3 0.3];
bin_size = 0.005;
convert_to_sec = 1;
sp_times_bin_norm = cell(numel(sp_times_filt), 1);
for i = 1:numel(sp_times_filt)
    % align and bin firing rates
    sp_times_bin_temp = get_bin_firing_rates(sp_times_filt{i}, sp_contact_filt{i}, lick_num, align_lick_num, smooth_bin_window, bin_size, convert_to_sec);
    
    % normalize firing rates to standard deviation as in Ames et al. 2019 
    std_temp = std(sp_times_bin_temp, [], 'all');
    if std_temp >= 1
        sp_times_bin_norm{i} = sp_times_bin_temp/std_temp;
    elseif std_temp < 1
        sp_times_bin_norm{i} = sp_times_bin_temp;
    end
end

% create PSTHs by averaging trials over conditions and smooth them
bin_size_back = 3;
bin_size_forward = 0;
mean_fakeout_firing_rate = nan(numel(sp_times_filt), 3, size(sp_times_bin_norm{1}, 2));
for i = 1:numel(sp_times_bin_norm)
    for j = 1:3
        mean_fakeout_firing_rate_temp = mean(sp_times_bin_norm{i}(fakeout_trial_filt{i} == j, :));
        mean_fakeout_firing_rate(i, j, :) = smoothdata(mean_fakeout_firing_rate_temp, 'movmean', [bin_size_back bin_size_forward]);
    end
end

% remove neurons whose firing rates equal zero
mean_fakeout_firing_rate = mean_fakeout_firing_rate(~mean(mean_fakeout_firing_rate, [2 3]) == 0, :, :);

% mean center each neuron's firing rate
mean_fakeout_firing_rate = mean_fakeout_firing_rate - mean(mean_fakeout_firing_rate, [2 3]);

% reshape the matrix for PCA, filter the bins that were for smooothing, then before PCA
bin_diff = ceil(abs(smooth_bin_window - pca_bin_window)/bin_size);
start_ind = 1 + bin_diff(1);
end_ind = size(mean_fakeout_firing_rate, 3) - bin_diff(2);
mean_fakeout_firing_rate_pca = [squeeze(mean_fakeout_firing_rate(:, 1, start_ind:end_ind)) squeeze(mean_fakeout_firing_rate(:, 2, start_ind:end_ind)) squeeze(mean_fakeout_firing_rate(:, 3, start_ind:end_ind))];
[coeff, ~, ~, ~, explained, ~] = pca(mean_fakeout_firing_rate_pca', 'Algorithm', 'svd', 'Centered', false, 'Economy', true);

% project the data onto top PCs
top_pcs = coeff(:, 1:3);
fakeout_pca = mean_fakeout_firing_rate_pca'*top_pcs;
% create the 3D plot, with dots indicated beginning of trajectories, L2
% contact, and L3 contact.
start_ind = 1;
end_ind = size(mean_fakeout_firing_rate_pca, 2)/3;
contact_bin = abs(pca_bin_window/bin_size) + 1;
L3_prot_time = round(mean(cellfun(@(x) mean(x(2, :)), prot_filt)));
L3_contact_time = round(mean(cellfun(@(x) mean(x(2, :)), sp_contact_filt)));
L2_contact_time = round(mean(cellfun(@(x) mean(x(1, :)), sp_contact_filt)));
L3_contact_bin = round((L3_contact_time - L2_contact_time)/(bin_size*1000));
L3_prot_bin = round((L3_prot_time - L2_contact_time)/(bin_size*1000));
h = figure();
for i = 1:size(mean_fakeout_firing_rate, 2)
    plot3(fakeout_pca(start_ind:end_ind, 1), fakeout_pca(start_ind:end_ind, 2), fakeout_pca(start_ind:end_ind, 3), 'linewidth', 5);
    hold on; grid on;
    scatter3(fakeout_pca(start_ind, 1), fakeout_pca(start_ind, 2), fakeout_pca(start_ind, 3), 100, 'filled', 'black');
    scatter3(fakeout_pca(contact_bin(1), 1), fakeout_pca(contact_bin(1), 2), fakeout_pca(contact_bin(1), 3), 100, 'filled', 'red');
    scatter3(fakeout_pca(L3_contact_bin + contact_bin(1), 1), fakeout_pca(L3_contact_bin + contact_bin(1), 2), fakeout_pca(L3_contact_bin + contact_bin(1), 3), 100, 'filled', 'blue');
    scatter3(fakeout_pca(L3_prot_bin + contact_bin(1), 1), fakeout_pca(L3_prot_bin + contact_bin(1), 2), fakeout_pca(L3_prot_bin + contact_bin(1), 3), 100, 'filled', 'yellow');
    start_ind = start_ind + size(mean_fakeout_firing_rate_pca, 2)/3;
    end_ind = end_ind + size(mean_fakeout_firing_rate_pca, 2)/3;
    contact_bin = contact_bin + size(mean_fakeout_firing_rate_pca, 2)/3;
end

xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
xlim([-7 7]);
ylim([-4 4]);
zlim([-5 5]);
xticks([-7 0 7]);
yticks([-4 0 4]);
zticks([-5 0 5]);

% Order of conditions is: 1) left, 2) center, 3) right 
dist_length = size(fakeout_pca, 1)/3;

left_center_dist = vecnorm([fakeout_pca(1:dist_length, :) - fakeout_pca((dist_length+1):(dist_length*2), :)]');
right_center_dist = vecnorm([fakeout_pca(((dist_length*2)+1):(dist_length*3), :) - fakeout_pca((dist_length+1):(dist_length*2), :)]');

left_center_dist_rand = nan(num_shuffles, numel(left_center_dist));
right_center_dist_rand = nan(num_shuffles, numel(right_center_dist));
z_left_center_dist_thresh = nan(num_shuffles, 1);
z_right_center_dist_thresh = nan(num_shuffles, 1);
for x = 1:num_shuffles
    % randomly select neurons with replacement
    neuron_ind = randi([1 numel(sp_times_bin_norm)], numel(sp_times_bin_norm), 1);
    
    % filter spike times and fakeout trials with these randomized indices
    sp_times_bin_norm_rand = sp_times_bin_norm(neuron_ind); 
    fakeout_trial_filt_rand = fakeout_trial_filt(neuron_ind);
    
    % randomly select trials for each neuron with replacement
    for y = 1:numel(sp_times_bin_norm_rand)
        trial_ind = randi([1 size(sp_times_bin_norm_rand{y}, 1)], size(sp_times_bin_norm_rand{y}, 1), 1);
        
        % filter spike times and fakeout trials with these randomized inds
        sp_times_bin_norm_rand{y} = sp_times_bin_norm_rand{y}(trial_ind, :);
        fakeout_trial_filt_rand{y} = fakeout_trial_filt_rand{y}(trial_ind);
    end
    
    % create PSTHs by averaging trials over conditions and smooth them
    bin_size_back = 3;
    bin_size_forward = 0;
    mean_fakeout_firing_rate_rand = nan(numel(sp_times_bin_norm_rand), 3, size(sp_times_bin_norm_rand{1}, 2));
    for i = 1:numel(sp_times_bin_norm_rand)
        for j = 1:3
            mean_fakeout_firing_rate_temp = mean(sp_times_bin_norm_rand{i}(fakeout_trial_filt_rand{i} == j, :));
            mean_fakeout_firing_rate_rand(i, j, :) = smoothdata(mean_fakeout_firing_rate_temp, 'movmean', [bin_size_back bin_size_forward]);
        end
    end

    % remove neurons whose firing rates equal zero
    mean_fakeout_firing_rate_rand = mean_fakeout_firing_rate_rand(~mean(mean_fakeout_firing_rate_rand, [2 3]) == 0, :, :);

    % mean center each neuron's firing rate
    mean_fakeout_firing_rate_rand = mean_fakeout_firing_rate_rand - mean(mean_fakeout_firing_rate_rand, [2 3]);

    % reshape the matrix for PCA, filter the bins that were for smooothing, then before PCA
    bin_diff = ceil(abs(smooth_bin_window - pca_bin_window)/bin_size);
    start_ind = 1 + bin_diff(1);
    end_ind = size(mean_fakeout_firing_rate_rand, 3) - bin_diff(2);
    mean_fakeout_firing_rate_pca_rand = [squeeze(mean_fakeout_firing_rate_rand(:, 1, start_ind:end_ind)) squeeze(mean_fakeout_firing_rate_rand(:, 2, start_ind:end_ind)) squeeze(mean_fakeout_firing_rate_rand(:, 3, start_ind:end_ind))];
    [coeff_rand, ~, ~, ~, ~, ~] = pca(mean_fakeout_firing_rate_pca_rand', 'Algorithm', 'svd', 'Centered', false, 'Economy', true);

    % project the data onto top PCs
    top_pcs_rand = coeff_rand(:, 1:19);
    fakeout_pca_rand = mean_fakeout_firing_rate_pca_rand'*top_pcs_rand;

    dist_length_temp = size(fakeout_pca_rand, 1)/3;
    left_center_dist_rand(x, :) = vecnorm([fakeout_pca_rand(1:dist_length_temp, :) - fakeout_pca_rand((dist_length_temp+1):(dist_length_temp*2), :)]');
    right_center_dist_rand(x, :) = vecnorm([fakeout_pca_rand(((dist_length_temp*2)+1):(dist_length_temp*3), :) - fakeout_pca_rand((dist_length_temp+1):(dist_length_temp*2), :)]');
    baseline_ind = size(left_center_dist_rand, 2)/2;
    z_left_center_dist_rand = (left_center_dist_rand(x, :) - mean(left_center_dist_rand(x, 1:baseline_ind)))/std(left_center_dist_rand(x, 1:baseline_ind));
    z_right_center_dist_rand = (right_center_dist_rand(x, :) - mean(right_center_dist_rand(x, 1:baseline_ind)))/std(right_center_dist_rand(x, 1:baseline_ind));
 
    z_left_center_dist_temp = find(z_left_center_dist_rand((baseline_ind+1):(baseline_ind*2)) > 3, 1);
    if ~isempty(z_left_center_dist_temp)
        z_left_center_dist_thresh(x) = z_left_center_dist_temp;
    else
        z_left_center_dist_thresh(x) = NaN;
    end

    z_right_center_dist_temp = find(z_right_center_dist_rand((baseline_ind+1):(baseline_ind*2)) > 3, 1);
    if ~isempty(z_right_center_dist_temp)
        z_right_center_dist_thresh(x) = z_right_center_dist_temp;
    else
        z_right_center_dist_thresh(x) = NaN;
    end
    
end

% Get median divergence time for each condition
med_left_center_div = median(z_left_center_dist_thresh, "omitnan"); %"omitmissing");
med_right_center_div = median(z_right_center_dist_thresh, "omitnan"); %"omitmissing");
med_center_div = median([z_left_center_dist_thresh; z_right_center_dist_thresh], "omitnan"); %"omitmissing");
iqr_center_div = prctile([z_left_center_dist_thresh; z_right_center_dist_thresh], [25 75]);

% Calculate SEM
std_left_center_div = std(left_center_dist_rand, 1);
std_right_center_div = std(right_center_dist_rand, 1);

figure();

% create SEM error shading for left distance
x_axis = linspace(1, size(left_center_dist_rand, 2), numel(left_center_dist));
std_upper = left_center_dist + std_left_center_div;
std_lower = left_center_dist - std_left_center_div;
x_axis2 = [x_axis, fliplr(x_axis)];
in_between = [std_upper, fliplr(std_lower)];
h = fill(x_axis2, in_between, 'r');
h.FaceAlpha = 0.2;
h.EdgeColor = 'none';
hold on;

% create SEM error shading for right distance
x_axis = linspace(1, size(left_center_dist_rand, 2), numel(right_center_dist));
std_upper = right_center_dist + std_right_center_div;
std_lower = right_center_dist - std_right_center_div;
x_axis2 = [x_axis, fliplr(x_axis)];
in_between = [std_upper, fliplr(std_lower)];
h = fill(x_axis2, in_between, 'c');
h.FaceAlpha = 0.2;
h.EdgeColor = 'none';

plot(left_center_dist); hold on;
plot(right_center_dist);
ylim([0 5]);

% plot divergence line - note that 61 is first timestep
med_div_both = length(left_center_dist)/2 + median([med_left_center_div med_right_center_div]);
plot([med_div_both med_div_both], [0 7]);

% plot L3 protrusion time - note that we are referencing to L2 contact time
L3_prot_temp = length(left_center_dist)/2 + L3_prot_bin;
plot([L3_prot_temp L3_prot_temp], [0 7]);

end