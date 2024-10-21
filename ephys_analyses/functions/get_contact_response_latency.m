function [wrs_hist, hist_bins, min_wrs_latency, min_wrs_latency_sig, min_wrs_latency_boot] = get_contact_response_latency(sp_vars_struct, selectivity, params, lick_num)

sp_times_filt = {sp_vars_struct.sp_times_filt};
sp_contact_filt = {sp_vars_struct.sp_contact_filt};
fakeout_trial_filt = {sp_vars_struct.fakeout_trial_filt};

[wrs_hist, hist_bins, min_wrs_latency, min_wrs_latency_sig] = calculate_contact_response_latency(sp_times_filt, sp_contact_filt, fakeout_trial_filt, selectivity, params, lick_num);

% Do a hierarchical bootstrap
min_wrs_latency_boot = nan(numel(sp_vars_struct), 1);
if params.bootstrap == true

    num_shuffles = params.num_shuffles;
    if params.par_shuf == true
        parfor i = 1:num_shuffles

            disp(i)
        
            rand_neuron_ind = randi(numel(sp_vars_struct), 1, numel(sp_vars_struct));
            sp_vars_struct_rand = sp_vars_struct(1, rand_neuron_ind);
            selectivity_rand = selectivity(rand_neuron_ind, :);
    
            sp_times_rand_temp = {sp_vars_struct_rand.sp_times_filt};
            sp_contact_rand_temp = {sp_vars_struct_rand.sp_contact_filt};
            fakeout_trial_rand_temp = {sp_vars_struct_rand.fakeout_trial_filt};
        
            sp_times_rand = cell(1, numel(sp_vars_struct));
            sp_contact_rand = cell(1, numel(sp_vars_struct));
            fakeout_trial_rand = cell(1, numel(sp_vars_struct));
            for j = 1:numel(sp_vars_struct_rand)
                num_trials_temp = numel(sp_vars_struct_rand(j).sp_times_filt);
                rand_trial_ind = randi(num_trials_temp, 1, num_trials_temp);
    
                sp_times_rand{1, j} = sp_times_rand_temp{j}(:, rand_trial_ind);
                sp_contact_rand{1, j} = sp_contact_rand_temp{j}(:, rand_trial_ind);
                fakeout_trial_rand{1, j} = fakeout_trial_rand_temp{j}(:, rand_trial_ind);
            end
    
        [~, ~, ~, min_wrs_latency_sig_temp] = calculate_contact_response_latency(sp_times_rand, sp_contact_rand, fakeout_trial_rand, selectivity_rand, params, lick_num);
        
        min_wrs_latency_boot(i, 1) = median(min_wrs_latency_sig_temp, "omitnan");
        end
    else
        for i = 1:num_shuffles

            disp(i)
        
            rand_neuron_ind = randi(numel(sp_vars_struct), 1, numel(sp_vars_struct));
            sp_vars_struct_rand = sp_vars_struct(1, rand_neuron_ind);
            selectivity_rand = selectivity(rand_neuron_ind, :);
    
            sp_times_rand_temp = {sp_vars_struct_rand.sp_times_filt};
            sp_contact_rand_temp = {sp_vars_struct_rand.sp_contact_filt};
            fakeout_trial_rand_temp = {sp_vars_struct_rand.fakeout_trial_filt};
        
            sp_times_rand = cell(1, numel(sp_vars_struct));
            sp_contact_rand = cell(1, numel(sp_vars_struct));
            fakeout_trial_rand = cell(1, numel(sp_vars_struct));
            for j = 1:numel(sp_vars_struct_rand)
                num_trials_temp = numel(sp_vars_struct_rand(j).sp_times_filt);
                rand_trial_ind = randi(num_trials_temp, 1, num_trials_temp);
    
                sp_times_rand{1, j} = sp_times_rand_temp{j}(:, rand_trial_ind);
                sp_contact_rand{1, j} = sp_contact_rand_temp{j}(:, rand_trial_ind);
                fakeout_trial_rand{1, j} = fakeout_trial_rand_temp{j}(:, rand_trial_ind);
            end
    
        [~, ~, ~, min_wrs_latency_sig_temp] = calculate_contact_response_latency(sp_times_rand, sp_contact_rand, fakeout_trial_rand, selectivity_rand, params, lick_num);
        
        min_wrs_latency_boot(1, i) = median(min_wrs_latency_sig_temp, "omitnan");
        end
    end
else
    min_wrs_latency_boot = NaN;
end


%% Function List

function [wrs_hist, hist_bins, min_wrs_latency, min_wrs_latency_sig] = calculate_contact_response_latency(sp_times_filt, sp_contact_filt, fakeout_trial_filt, selectivity, params, lick_num)

% note that, we want to have spare bins to the left and right of what we
% care about in bin_window for smoothing of PSTHs below.
align_lick_num = params.align_lick_num;
bin_window = params.bin_window;
bin_size = params.bin_size;
convert_to_sec = 1;
sp_times_bin_all = cell(numel(sp_times_filt), 1);
sp_times_bin_cond = cell(size(sp_times_bin_all, 1), 3);
for i = 1:numel(sp_times_filt)
    % align and bin firing rates
    sp_times_bin_temp = get_bin_firing_rates(sp_times_filt{i}, sp_contact_filt{i}, lick_num, align_lick_num, bin_window, bin_size, convert_to_sec);
    sp_times_bin_all{i} = sp_times_bin_temp;
    for j = 1:3
        fakeout_firing_rate_temp = sp_times_bin_temp(fakeout_trial_filt{i} == j, :);
        sp_times_bin_cond{i, j} = fakeout_firing_rate_temp;
    end
end

% number of bins to take sum of spikes in - including current bin. e.g.,
% the number of bins is + 1 of the entry.
bin_size_back = params.bin_size_back;
bin_size_forward = params.bin_size_forward;
bin_step = params.bin_step;
start_bin = size(sp_times_bin_all{1}, 2)/2 + 1 + bin_size_back;
bin_ind = start_bin:bin_step:size(sp_times_bin_all{1}, 2);
wrs_pvals = nan(size(sp_times_bin_cond, 1), size(sp_times_bin_cond, 2), numel(bin_ind));
wrs_stats = cell(size(sp_times_bin_cond, 1), size(sp_times_bin_cond, 2), numel(bin_ind));
for i = 1:size(sp_times_bin_cond, 1)
    k = 1;
    for j = bin_ind
        sum_spikes_1 = sum(sp_times_bin_cond{i, 1}(:, (j - bin_size_back):(j + bin_size_forward)), 2);
        sum_spikes_2 = sum(sp_times_bin_cond{i, 2}(:, (j - bin_size_back):(j + bin_size_forward)), 2);
        sum_spikes_3 = sum(sp_times_bin_cond{i, 3}(:, (j - bin_size_back):(j + bin_size_forward)), 2);

        [wrs_pvals(i, 1, k), ~, wrs_stats{i, 1, k}] = ranksum(sum_spikes_1, sum_spikes_2);
        [wrs_pvals(i, 2, k), ~, wrs_stats{i, 2, k}] = ranksum(sum_spikes_3, sum_spikes_2);
        [wrs_pvals(i, 3, k), ~, wrs_stats{i, 3, k}] = ranksum(sum_spikes_1, sum_spikes_3);
        k = k + 1;
    end
end

num_bins_sig = params.num_bins_sig;
pval = params.pval;
wrs_pvals_bin = wrs_pvals < pval;

wrs_latency = nan(size(wrs_pvals, 1), size(wrs_pvals, 2));
for i = 1:size(wrs_pvals, 1)
    for j = 1:size(wrs_pvals, 2)
        temp = squeeze(wrs_pvals_bin(i, j, :));
        temp_ind = find(temp);
        temp_ind_diff = diff(temp_ind) == 1; % diff of one indicates consecutive indices
        num_bins_sig_ind = num_bins_sig - 1;
        for k = 1:numel(temp_ind_diff) - num_bins_sig_ind
            sum_sig_bins = sum(temp_ind_diff(k:k+num_bins_sig_ind));
            if sum_sig_bins == num_bins_sig
                wrs_latency(i, j) = temp_ind(k)*bin_step;
                break
            end
        end
    end
end

min_wrs_latency = min(wrs_latency, [], 2);
selectivity_L2 = [selectivity{:, 1}]';
min_wrs_latency_sig = min_wrs_latency(selectivity_L2 ~= 0);

ignore_first_bin = params.ignore_first_bin;
if ignore_first_bin == true
    min_wrs_latency_sig = min_wrs_latency_sig(min_wrs_latency_sig > (bin_size_back + 1));

end

hist_bin_size = bin_size_back+1;
hist_bins = 1:hist_bin_size:(bin_window(2)*1000);
wrs_hist = nan(numel(hist_bins) - 1, 1);
for i = 1:numel(hist_bins) - 1
    wrs_hist(i, 1) = sum(min_wrs_latency_sig >= hist_bins(i) & min_wrs_latency_sig < hist_bins(i+1));

end