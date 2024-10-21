function create_selectivity_heatmaps(sp_vars_struct, sig_struct, params, lick_num)

sp_times_filt = {sp_vars_struct.sp_times_filt};
sp_contact_filt = {sp_vars_struct.sp_contact_filt};
fakeout_trial_filt = {sp_vars_struct.fakeout_trial_filt};

align_lick_num = params.align_lick_num;
smooth_bin_window = params.smooth_bin_window;
real_window = params.real_window;
bin_size = params.bin_size;
bin_size_back = params.bin_size_back;
bin_size_forward = params.bin_size_forward;

% note that, we want to have spare bins to the left and right of what we
% care about in bin_window for smoothing of PSTHs below.
filt_window = (smooth_bin_window - real_window)/bin_size;
filt_ind = round([abs(filt_window(1)) + 1; sum(abs(smooth_bin_window))/bin_size - abs(filt_window(1))]);
convert_to_sec = 1;
sp_times_bin_norm = cell(numel(sp_vars_struct), 1);
for i = 1:numel(sp_vars_struct)
    % align and bin firing rates
    sp_times_bin_temp = get_bin_firing_rates(sp_times_filt{i}, sp_contact_filt{i}, lick_num, align_lick_num, smooth_bin_window, bin_size, convert_to_sec);
    sp_times_bin_norm{i}  = sp_times_bin_temp;
end

% create PSTHs by averaging trials over conditions and smooth them
mean_fakeout_firing_rate = nan(numel(sp_vars_struct), 3, size(sp_times_bin_norm{1}, 2) - sum(abs(filt_window)));
for i = 1:numel(sp_times_bin_norm)
    for j = 1:3
        mean_fakeout_firing_rate_temp = mean(sp_times_bin_norm{i}(fakeout_trial_filt{i} == j, :));
        mean_fakeout_firing_rate_temp = smoothdata(mean_fakeout_firing_rate_temp, 'movmean', [bin_size_back bin_size_forward]);
        mean_fakeout_firing_rate(i, j, :) = mean_fakeout_firing_rate_temp(filt_ind(1):filt_ind(2));
    end
end

% determine which neurons have contra/ipsi/center/etc selectivity
selectivity_cent = sig_struct.selectivity_cent;
session_name = {selectivity_cent{:, 3}};
sig_left_contra = cell2mat(cellfun(@(x) contains(x, 'right'), session_name, 'UniformOutput', false)) & [selectivity_cent{:, 1}] == 1;
sig_right_contra = cell2mat(cellfun(@(x) contains(x, 'left'), session_name, 'UniformOutput', false)) & [selectivity_cent{:, 1}] == 3;
sig_contra = sig_left_contra | sig_right_contra;

sig_right_ipsi = cell2mat(cellfun(@(x) contains(x, 'right'), session_name, 'UniformOutput', false)) & [selectivity_cent{:, 1}] == 3;
sig_left_ipsi = cell2mat(cellfun(@(x) contains(x, 'left'), session_name, 'UniformOutput', false)) & [selectivity_cent{:, 1}] == 1;
sig_ipsi = sig_left_ipsi | sig_right_ipsi;

sig_right_center= cell2mat(cellfun(@(x) contains(x, 'right'), session_name, 'UniformOutput', false)) & [selectivity_cent{:, 1}] == 2;
sig_left_center = cell2mat(cellfun(@(x) contains(x, 'left'), session_name, 'UniformOutput', false)) & [selectivity_cent{:, 1}] == 2;
sig_center = sig_left_center | sig_right_center;

sig_right_contra_contracenter= cell2mat(cellfun(@(x) contains(x, 'right'), session_name, 'UniformOutput', false)) & [selectivity_cent{:, 1}] == 4;
sig_left_contra_contracenter = cell2mat(cellfun(@(x) contains(x, 'left'), session_name, 'UniformOutput', false)) & [selectivity_cent{:, 1}] == 5;
sig_contra_and_center = sig_left_contra_contracenter | sig_right_contra_contracenter;

sig_right_ipsi_ipsicenter= cell2mat(cellfun(@(x) contains(x, 'right'), session_name, 'UniformOutput', false)) & [selectivity_cent{:, 1}] == 5;
sig_left_ipsi_ipsicenter = cell2mat(cellfun(@(x) contains(x, 'left'), session_name, 'UniformOutput', false)) & [selectivity_cent{:, 1}] == 4;
sig_ipsi_and_center = sig_left_ipsi_ipsicenter | sig_right_ipsi_ipsicenter;

% min-max normalization of the firing rates
mean_fakeout_firing_rate_filt = nan(numel(sp_vars_struct), 3, size(sp_times_bin_norm{1}, 2) - sum(abs(filt_window)));
for i = 1:size(mean_fakeout_firing_rate, 1)
    min_fr = min(mean_fakeout_firing_rate(i, :, :), [], 'all');
    max_fr = max(mean_fakeout_firing_rate(i, :, :), [], 'all');
    mean_fakeout_firing_rate_filt(i, :, :) = (mean_fakeout_firing_rate(i, :, :) - min_fr)/(max_fr - min_fr);
end

% Construct the matrices for heatmaps

% Contra Selective Neurons
k = 1;
mean_fr_contra = nan(sum(sig_contra), size(mean_fakeout_firing_rate, 3));
mean_fr_center = nan(sum(sig_contra), size(mean_fakeout_firing_rate, 3));
mean_fr_ipsi = nan(sum(sig_contra), size(mean_fakeout_firing_rate, 3));
max_fr_contra_ind = nan(sum(sig_contra), 1);
align_time_ind = size(mean_fakeout_firing_rate, 3)/2 + 1;
for i = 1:numel(sig_left_contra)
    if sig_left_contra(i) == 1
        mean_fr_contra(k, :) = squeeze(mean_fakeout_firing_rate_filt(i, 1, :)); 
        mean_fr_center(k, :) = squeeze(mean_fakeout_firing_rate_filt(i, 2, :)); 
        mean_fr_ipsi(k, :) = squeeze(mean_fakeout_firing_rate_filt(i, 3, :)); 
        max_fr_temp = find(mean_fr_contra(k, align_time_ind:end) == max(mean_fr_contra(k, align_time_ind:end)));
        max_fr_contra_ind(k, 1) = max_fr_temp(1);
        k = k + 1;
    elseif sig_right_contra(i) == 1
        mean_fr_contra(k, :) = squeeze(mean_fakeout_firing_rate_filt(i, 3, :)); 
        mean_fr_center(k, :) = squeeze(mean_fakeout_firing_rate_filt(i, 2, :)); 
        mean_fr_ipsi(k, :) = squeeze(mean_fakeout_firing_rate_filt(i, 1, :)); 
        max_fr_temp = find(mean_fr_contra(k, align_time_ind:end) == max(mean_fr_contra(k, align_time_ind:end)));
        max_fr_contra_ind(k, 1) = max_fr_temp(1);
        k = k + 1;
    end
end

[~, max_fr_ind] = sort(max_fr_contra_ind);

figure(); imagesc([mean_fr_contra(max_fr_ind, :); mean_fr_center(max_fr_ind, :); mean_fr_ipsi(max_fr_ind, :)]);
yticks([sum(sig_contra), sum(sig_contra)*2, sum(sig_contra)*3])
xticks([1 align_time_ind size(mean_fakeout_firing_rate, 3)])
xticklabels({'-300', '0', '300'})
xlim([1 size(mean_fakeout_firing_rate, 3)])
colorbar
clim([0 1])
title('Contra Selective')

% Ipsi Selective Neurons
k = 1;
mean_fr_contra = nan(sum(sig_ipsi), size(mean_fakeout_firing_rate, 3));
mean_fr_center = nan(sum(sig_ipsi), size(mean_fakeout_firing_rate, 3));
mean_fr_ipsi = nan(sum(sig_ipsi), size(mean_fakeout_firing_rate, 3));
max_fr_ipsi_ind = nan(sum(sig_ipsi), 1);
align_time_ind = size(mean_fakeout_firing_rate, 3)/2 + 1;
for i = 1:numel(sig_left_ipsi)
    if sig_left_ipsi(i) == 1
        mean_fr_contra(k, :) = squeeze(mean_fakeout_firing_rate_filt(i, 3, :)); 
        mean_fr_center(k, :) = squeeze(mean_fakeout_firing_rate_filt(i, 2, :)); 
        mean_fr_ipsi(k, :) = squeeze(mean_fakeout_firing_rate_filt(i, 1, :)); 
        max_fr_temp = find(mean_fr_ipsi(k, align_time_ind:end) == max(mean_fr_ipsi(k, align_time_ind:end)));
        max_fr_ipsi_ind(k, 1) = max_fr_temp(1);
        k = k + 1;
    elseif sig_right_ipsi(i) == 1
        mean_fr_contra(k, :) = squeeze(mean_fakeout_firing_rate_filt(i, 1, :)); 
        mean_fr_center(k, :) = squeeze(mean_fakeout_firing_rate_filt(i, 2, :)); 
        mean_fr_ipsi(k, :) = squeeze(mean_fakeout_firing_rate_filt(i, 3, :)); 
        max_fr_temp = find(mean_fr_ipsi(k, align_time_ind:end) == max(mean_fr_ipsi(k, align_time_ind:end)));
        max_fr_ipsi_ind(k, 1) = max_fr_temp(1);
        k = k + 1;
    end
end

[~, max_fr_ind] = sort(max_fr_ipsi_ind);

figure(); imagesc([mean_fr_contra(max_fr_ind, :); mean_fr_center(max_fr_ind, :); mean_fr_ipsi(max_fr_ind, :)]);
yticks([sum(sig_ipsi), sum(sig_ipsi)*2, sum(sig_ipsi)*3])
xticks([1 align_time_ind size(mean_fakeout_firing_rate, 3)])
xticklabels({'-300', '0', '300'})
xlim([1 size(mean_fakeout_firing_rate, 3)])
colorbar
clim([0 1])
title('Ipsi Selective')

% Center Selective Neurons
k = 1;
mean_fr_contra = nan(sum(sig_center), size(mean_fakeout_firing_rate, 3));
mean_fr_center = nan(sum(sig_center), size(mean_fakeout_firing_rate, 3));
mean_fr_ipsi = nan(sum(sig_center), size(mean_fakeout_firing_rate, 3));
max_fr_center_ind = nan(sum(sig_center), 1);
align_time_ind = size(mean_fakeout_firing_rate, 3)/2 + 1;
for i = 1:numel(sig_left_center)
    if sig_left_center(i) == 1
        mean_fr_contra(k, :) = squeeze(mean_fakeout_firing_rate_filt(i, 3, :)); 
        mean_fr_center(k, :) = squeeze(mean_fakeout_firing_rate_filt(i, 2, :)); 
        mean_fr_ipsi(k, :) = squeeze(mean_fakeout_firing_rate_filt(i, 1, :)); 
        max_fr_temp = find(mean_fr_center(k, align_time_ind:end) == max(mean_fr_center(k, align_time_ind:end)));
        max_fr_center_ind(k, 1) = max_fr_temp(1);
        k = k + 1;
    elseif sig_right_center(i) == 1
        mean_fr_contra(k, :) = squeeze(mean_fakeout_firing_rate_filt(i, 1, :)); 
        mean_fr_center(k, :) = squeeze(mean_fakeout_firing_rate_filt(i, 2, :)); 
        mean_fr_ipsi(k, :) = squeeze(mean_fakeout_firing_rate_filt(i, 3, :)); 
        max_fr_temp = find(mean_fr_center(k, align_time_ind:end) == max(mean_fr_center(k, align_time_ind:end)));
        max_fr_center_ind(k, 1) = max_fr_temp(1);
        k = k + 1;
    end
end

[~, max_fr_ind] = sort(max_fr_center_ind);

figure(); imagesc([mean_fr_contra(max_fr_ind, :); mean_fr_center(max_fr_ind, :); mean_fr_ipsi(max_fr_ind, :)]);
yticks([sum(sig_center), sum(sig_center)*2, sum(sig_center)*3])
xticks([1 align_time_ind size(mean_fakeout_firing_rate, 3)])
xticklabels({'-300', '0', '300'})
xlim([1 size(mean_fakeout_firing_rate, 3)])
colorbar
clim([0 1])
title('Center Selective')

% Contra + Center Selective Neurons
k = 1;
mean_fr_contra = nan(sum(sig_contra_and_center), size(mean_fakeout_firing_rate, 3));
mean_fr_center = nan(sum(sig_contra_and_center), size(mean_fakeout_firing_rate, 3));
mean_fr_ipsi = nan(sum(sig_contra_and_center), size(mean_fakeout_firing_rate, 3));
max_fr_center_ind = nan(sum(sig_contra_and_center), 1);
align_time_ind = size(mean_fakeout_firing_rate, 3)/2 + 1;
for i = 1:numel(sig_left_contra_contracenter)
    if sig_left_contra_contracenter(i) == 1
        mean_fr_contra(k, :) = squeeze(mean_fakeout_firing_rate_filt(i, 3, :)); 
        mean_fr_center(k, :) = squeeze(mean_fakeout_firing_rate_filt(i, 2, :)); 
        mean_fr_ipsi(k, :) = squeeze(mean_fakeout_firing_rate_filt(i, 1, :)); 
        max_fr_temp = find(mean_fr_contra(k, align_time_ind:end) == max(mean_fr_contra(k, align_time_ind:end)));
        max_fr_center_ind(k, 1) = max_fr_temp(1);
        k = k + 1;
    elseif sig_right_contra_contracenter(i) == 1
        mean_fr_contra(k, :) = squeeze(mean_fakeout_firing_rate_filt(i, 1, :)); 
        mean_fr_center(k, :) = squeeze(mean_fakeout_firing_rate_filt(i, 2, :)); 
        mean_fr_ipsi(k, :) = squeeze(mean_fakeout_firing_rate_filt(i, 3, :)); 
        max_fr_temp = find(mean_fr_contra(k, align_time_ind:end) == max(mean_fr_contra(k, align_time_ind:end)));
        max_fr_center_ind(k, 1) = max_fr_temp(1);
        k = k + 1;
    end
end

[~, max_fr_ind] = sort(max_fr_center_ind);

figure(); imagesc([mean_fr_contra(max_fr_ind, :); mean_fr_center(max_fr_ind, :); mean_fr_ipsi(max_fr_ind, :)]);
yticks([sum(sig_contra_and_center), sum(sig_contra_and_center)*2, sum(sig_contra_and_center)*3])
xticks([1 align_time_ind size(mean_fakeout_firing_rate, 3)])
xticklabels({'-300', '0', '300'})
xlim([1 size(mean_fakeout_firing_rate, 3)])
colorbar
clim([0 1])
title('Contra + Center Selective')

% Ipsi + Center Selective Neurons
sig_right_ipsi_ipsicenter= cell2mat(cellfun(@(x) contains(x, 'right'), session_name, 'UniformOutput', false)) & [selectivity_cent{:, 1}] == 5;
sig_left_ipsi_ipsicenter = cell2mat(cellfun(@(x) contains(x, 'left'), session_name, 'UniformOutput', false)) & [selectivity_cent{:, 1}] == 4;
sig_ipsi_and_center = sig_left_ipsi_ipsicenter | sig_right_ipsi_ipsicenter;

k = 1;
mean_fr_contra = nan(sum(sig_ipsi_and_center), size(mean_fakeout_firing_rate, 3));
mean_fr_center = nan(sum(sig_ipsi_and_center), size(mean_fakeout_firing_rate, 3));
mean_fr_ipsi = nan(sum(sig_ipsi_and_center), size(mean_fakeout_firing_rate, 3));
max_fr_ipsi_ind = nan(sum(sig_ipsi_and_center), 1);
align_time_ind = size(mean_fakeout_firing_rate, 3)/2 + 1;
for i = 1:numel(sig_left_ipsi_ipsicenter)
    if sig_left_ipsi_ipsicenter(i) == 1
        mean_fr_contra(k, :) = squeeze(mean_fakeout_firing_rate_filt(i, 3, :)); 
        mean_fr_center(k, :) = squeeze(mean_fakeout_firing_rate_filt(i, 2, :)); 
        mean_fr_ipsi(k, :) = squeeze(mean_fakeout_firing_rate_filt(i, 1, :)); 
        max_fr_temp = find(mean_fr_ipsi(k, align_time_ind:end) == max(mean_fr_ipsi(k, align_time_ind:end)));
        max_fr_ipsi_ind(k, 1) = max_fr_temp(1);
        k = k + 1;
    elseif sig_right_ipsi_ipsicenter(i) == 1
        mean_fr_contra(k, :) = squeeze(mean_fakeout_firing_rate_filt(i, 1, :)); 
        mean_fr_center(k, :) = squeeze(mean_fakeout_firing_rate_filt(i, 2, :)); 
        mean_fr_ipsi(k, :) = squeeze(mean_fakeout_firing_rate_filt(i, 3, :)); 
        max_fr_temp = find(mean_fr_ipsi(k, align_time_ind:end) == max(mean_fr_ipsi(k, align_time_ind:end)));
        max_fr_ipsi_ind(k, 1) = max_fr_temp(1);
        k = k + 1;
    end
end

[~, max_fr_ind] = sort(max_fr_ipsi_ind);

figure(); imagesc([mean_fr_contra(max_fr_ind, :); mean_fr_center(max_fr_ind, :); mean_fr_ipsi(max_fr_ind, :)]);
yticks([sum(sig_ipsi_and_center), sum(sig_ipsi_and_center)*2, sum(sig_ipsi_and_center)*3])
xticks([1 align_time_ind size(mean_fakeout_firing_rate, 3)])
xticklabels({'-300', '0', '300'})
xlim([1 size(mean_fakeout_firing_rate, 3)])
colorbar
clim([0 1])
title('Center + Ipsi Selective')

end