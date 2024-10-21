%% Load full latSC dataset. All neurons are used for L2 analyses. 

load('latSC_ephys.mat')

%% Filter sp_struct to only include neurons recorded during re-centering sessions.

% Note this dataset is the only dataset where recording site was tracked.
session_id = {sp_struct.session_id}';
is_recenter = false(numel(sp_struct), 1);
for i = 1:numel(sp_struct)
    is_recenter_temp = regexp(session_id{i, 1}, '(?<=recentering).*');
    if isempty(is_recenter_temp)
        is_recenter(i, 1) = false;
    else
        is_recenter(i, 1) = true;
    end
end

sp_struct = sp_struct(is_recenter);

clearvars is_recenter_temp is_recenter i session_id

%% This packet of code is solely to reproduce the selectivity across AP axis within SC (Fig. 3c)

% Remove neurons with trials < min_trial_num & filter behavior vars

clearvars -except sp_struct

% replace with appropriate path to Accessory Files
load('Z:\bsi8\Ito, Gao et. al (2024)\Code Review\Accessory Files\site_map_recenter_SC.mat');

% this boolean indicates whether we are not (true) or are (false) analyzing
% re-centering sessions. In functions below, this is used to perform
% specific analyses/filtering. Yes, this should be changed to 'recenter' 
% instead of 'no_recenter', but I will change this later...
no_recenter = true;

% define lick numbers to be analyzed and params used for filtering trials.
% These analyses are focused on L2 and L3. We require L1 to L3 to exist
% (lick_exist_lick_num), to make spout contact on each lick
% (sp_contact_lick_num) only once (sp_contact2_absent_lick_num), and have
% contact tracking data (contact_centroid_lick_num).  
lick_num = [1 2 3];
params.lick_exist_lick_num = [1 2 3];
params.sp_contact_lick_num = [1 2 3];
params.sp_contact2_absent_lick_num = [1 2 3];
params.contact_centroid_lick_num = [1 2 3];

% exclude neurons with trials < min_trial_num & firing rate < min_fr_thresh,
% for trials that pass the behavioral params set above. 
min_trial_num = 10;
min_fr_thresh = 0.1;
[sp_vars_struct, sp_struct_filt] = filter_sp_neurons_and_vars_recenter(sp_struct, min_trial_num, min_fr_thresh, lick_num, params, no_recenter);

clearvars min_trial_num min_fr_thresh params

% Test for significance in fakeout and recenter conditions

% fakeout_trial_compare defines which licks to significance test across L2
% displacement conditions. spike_window defines the time window (sec) to
% analyze spiking activity in. fakeout_alpha defines Bonerroni-adj p-val.
params.fakeout_trial_compare = [1 3];
params.spike_window = [0 0.1];
num_comparisons = 3;
params.fakeout_alpha = 0.05/(num_comparisons);

% sig_lick_num defines the lick number to assess significance in firing
% rates. 
params.sig_lick_num = 2;
sig_struct = sig_test_firing_rates_fakeout_recenter(sp_struct_filt, sp_vars_struct, params, lick_num, no_recenter);

clearvars num_comparisons params

% Get recording site for selectivity - note that this only works with recentering dataset and with no_recenter == true

% assign recording site in latSC to selectivity
selectivity_site = assign_site_selectivity(sig_struct, site_map_recenter_SC);

num_shuffles = 100000;
[selectivity_bin, selectivity_rand_bin] = calculate_percent_site_selectivity(selectivity_site, num_shuffles);

site_sel_sig = plot_site_selectivity(selectivity_bin, selectivity_rand_bin, num_shuffles);

clearvars num_shuffles selectivity_bin selectivity_rand_bin site_map_recenter_SC selectivity_site

%% Now, we will filter the dataset for re-centering analyses. 

clearvars -except sp_struct

% this boolean indicates whether we are not (true) or are (false) analyzing
% re-centering sessions. In functions below, this is used to perform
% specific analyses/filtering. Yes, this should be changed to 'recenter' 
% instead of 'no_recenter', but I will change this later...
no_recenter = false;

% define lick numbers to be analyzed and params used for filtering trials.
% As these analyses are focused on L2 and L4. We require L1 to L5 to exist
% (lick_exist_lick_num), to make spout contact on each lick
% (sp_contact_lick_num) only once (sp_contact2_absent_lick_num), and have
% contact tracking data (contact_centroid_lick_num).  
lick_num = [1 2 3 4 5];
params.lick_exist_lick_num = [1 2 3 4 5];
params.sp_contact_lick_num = [1 2 3 4 5];
params.sp_contact2_absent_lick_num = [1 2 3 4 5];
params.contact_centroid_lick_num = [1 2 3 4 5];

% exclude neurons with trials < min_trial_num & firing rate < min_fr_thresh,
% for trials that pass the behavioral params set above. 
min_trial_num = 10;
min_fr_thresh = 0.1;
[sp_vars_struct, sp_struct_filt] = filter_sp_neurons_and_vars_recenter(sp_struct, min_trial_num, min_fr_thresh, lick_num, params, no_recenter);

clearvars min_trial_num min_fr_thresh params

%% Test for significance in fakeout and recenter conditions

% fakeout_trial_compare defines which licks to significance test across L2
% displacement conditions. spike_window defines the time window (sec) to
% analyze spiking activity in. fakeout_alpha defines Bonerroni-adj p-val.
params.fakeout_trial_compare = [1 3];
params.spike_window = [0 0.1];
num_comparisons = 3;
params.fakeout_alpha = 0.05/num_comparisons;

% sig_lick_num defines the lick number to assess significance in firing
% rates. We are also testing for re-centering significance.
params.sig_lick_num = [2 4];
params.recenter_alpha = 0.05/2;

sig_struct = sig_test_firing_rates_fakeout_recenter(sp_struct_filt, sp_vars_struct, params, lick_num, no_recenter);

clearvars num_comparisons

%% Plot percentage significant for fakeout, recenter, and fakeout and recenter. Extended Data Fig. 8c.

plot_percent_fakeout_recenter(sig_struct, params.sig_lick_num)

clearvars params

%% Plot percentage of fakeout/recenter neurons that were ACTIVATED by the same nick site

plot_same_nick_site_fakeout_recenter(sig_struct)

%% Calculate error vector and associated variables

% replace with appropriate path to Accessory Files
load('Z:\bsi8\Ito, Gao et. al (2024)\Code Review\Accessory Files\recenter_fiducial_list.mat')
contact_vars_struct = calculate_contact_error_vars(sp_vars_struct, fiducial_xy);

clearvars fiducial_xy

%% Calculate binned contact angle and firing rate correlations relative to tongue and head

% filter vars for only neurons modulated both in fakeout and recentering
fakeout_recenter_sig_ind = [selectivity{:, 1}] ~= 0 & ([selectivity{:, 2}] ~= 0 | [selectivity{:, 3}] ~= 0);

contact_ang_sig_hc = {contact_vars_struct(fakeout_recenter_sig_ind).contact_ang_hc};
contact_ang_sig = {contact_vars_struct(fakeout_recenter_sig_ind).contact_ang};
firing_rate_sig = sig_struct.firing_rate(fakeout_recenter_sig_ind, :, :);

% define the minimum and maximum angle values for the sliding window, as
% well as the bin width - create an array with start/end values
ang_min = -18;
ang_max = 18;
ang_bin_width = 3;
ang_bin = ang_min:ang_bin_width:ang_max;

% define minimum trial number needed to average
min_trial_num = 2;

% calculate the mean angles and mean firing rates 
num_shuffles = 1000;
[contact_angle_bin, firing_rate_bin, firing_rate_bin_shuff, num_trials_bin, mean_fr_shuff, mean_fr] = bin_behavior_and_firing_rate(contact_ang_sig, firing_rate_sig, ang_bin, min_trial_num, num_shuffles);
[contact_angle_bin_hc, firing_rate_bin_hc, firing_rate_bin_shuff_hc, num_trials_bin_hc, mean_fr_shuff_hc, mean_fr_hc] = bin_behavior_and_firing_rate(contact_ang_sig_hc, firing_rate_sig, ang_bin, min_trial_num, num_shuffles);

clearvars contact_ang_sig_hc contact_ang_sig min_trial_num ang_min ang_max ang_bin_width ang_bin x

%% Contact angle and firing rate correlation with shuffling for significance

% store p-values and corr. coefficients of each neuron. plot if desired
corr_lick_num = 4;
num_shuffles = 1000;
corr_r = nan(2, size(contact_angle_bin, 1));
corr_p = nan(2, size(contact_angle_bin, 1));
corr_coef_shuf = nan(2, size(contact_angle_bin, 1), num_shuffles);
p_shuff = nan(2, size(contact_angle_bin, 1));
r2_shuff_temp = nan(2, size(contact_angle_bin, 1), num_shuffles);

for x = 1:2
    if x == 1
        contact_angle_bin_temp = contact_angle_bin;
        firing_rate_bin_temp = firing_rate_bin;
        firing_rate_bin_shuff_temp = firing_rate_bin_shuff;
    elseif x == 2
        contact_angle_bin_temp = contact_angle_bin_hc;
        firing_rate_bin_temp = firing_rate_bin_hc;
        firing_rate_bin_shuff_temp = firing_rate_bin_shuff_hc;
    end

    for i = 1:size(contact_angle_bin_temp, 1)
    
        [r, p] = corrcoef(-13.5:3:13.5, firing_rate_bin_temp(i, lick_num == corr_lick_num, 2:11), 'rows', 'complete');
        corr_r(x, i) = r(1, 2).^2;
        corr_p(x, i) = p(1, 2);
    
        for j = 1:num_shuffles
            temp = corrcoef(-13.5:3:13.5, firing_rate_bin_shuff_temp(i, lick_num == corr_lick_num, 2:11, j), 'rows', 'complete');
            r2_shuff_temp(x, i, j) = temp(1, 2).^2;
            corr_coef_shuf(x, i, j) = temp(1, 2);
        end
        
        corr_95 = prctile(r2_shuff_temp(x, i, :), 95);
        p_shuff(x, i) = corr_r(x, i) > corr_95;
    
    end
end

% % plot a single neuron's correlation
% neuron = 1;
% figure('units','normalized','outerposition',[0 0 1 1]); 
% scatter(-13.5:3:13.5, squeeze(firing_rate_bin(neuron, lick_num == corr_lick_num, 2:11)));
% xlim([-15 15])

clearvars plot_neuron square_r x contact_angle_bin_temp firing_rate_bin_temp firing_rate_bin_shuff_temp i...
    r p j temp corr_95 corr_lick_num

%% Plot scatter of correlation coefficient tongue vs head-centered - Fig. 4e

tc = p_shuff(1, :) == 1 & corr_p(1, :) < 0.05 & (p_shuff(2, :) == 0 | corr_p(2, :) >= 0.05);
corr_r_tc = corr_r(:, tc); %30 36
hc = p_shuff(2, :) == 1 & corr_p(2, :) < 0.05 & (p_shuff(1, :) == 0 | corr_p(1, :) >= 0.05);
corr_r_hc = corr_r(:, hc); %21 14
both = p_shuff(1, :) == 1 & p_shuff(2, :) == 1 & corr_p(1, :) < 0.05 & corr_p(2, :) < 0.05;
corr_r_both = corr_r(:, both); %12 19
corr_r_not_sig = corr_r(:, (p_shuff(1, :) == 0 & p_shuff(2, :) == 0) | (corr_p(1, :) >= 0.05 & corr_p(2, :) >= 0.05)); %30 24

figure('units','normalized','outerposition',[0 0 1 1]);
scatter(corr_r_tc(1, :), corr_r_tc(2, :), 'filled', 'red'); hold on;
scatter(corr_r_hc(1, :), corr_r_hc(2, :), 'filled', 'black');
scatter(corr_r_both(1, :), corr_r_both(2, :), 'filled', 'blue');
scatter(corr_r_not_sig(1, :), corr_r_not_sig(2, :));
xlim([0 1]);
ylim([0 1]);

% Create filtered sp_vars for below
sp_vars_struct_sig = sp_vars_struct(1, fakeout_recenter_sig_ind);
selectivity_sig = selectivity(fakeout_recenter_sig_ind', 1);

sp_vars_struct_tc = sp_vars_struct_sig(1, tc);
selectivity_tc = selectivity_sig(tc, :);

sp_vars_struct_hc = sp_vars_struct_sig(1, hc);
selectivity_hc = selectivity_sig(hc, :);

sp_vars_struct_both = sp_vars_struct_sig(1, both);
selectivity_both = selectivity_sig(both, :);

sp_vars_struct_none = sp_vars_struct_sig(1, ~(tc | hc | both));
selectivity_none = selectivity_sig(~(tc | hc | both), :);

%% Plot percentage tongue, head, or tongue & head centric - Fig. 4e

hc = size(corr_r_hc, 2);
tc = size(corr_r_tc, 2);
hc_tc = size(corr_r_both, 2);
none = size(corr_r_not_sig, 2);
figure('units','normalized','outerposition',[0 0 1 1]);
mat_version = version;
if str2double(mat_version(1:2)) >= 23
    % MATLAB 2023B AND BEYOND: Bar chart of percentage contra vs. ispi fakeout neurons
    names = categorical({strcat('TC n=', num2str(tc)), strcat('HC n=', num2str(hc)), strcat('Both n=', num2str(hc_tc)), strcat('None n=', num2str(none))});
    donutchart([tc, hc, hc_tc, none], names);
else
    % PRE-MATLAB 2023B: Bar chart of percentage contra vs. ispi fakeout neurons
    names = [strcat("TC n=", num2str(tc)), strcat("HC n=", num2str(hc)), strcat("Both n=", num2str(hc_tc)), strcat("None n=", num2str(none))];
    pie([tc, hc, hc_tc, none], names);
end

clearvars hc tc hc_tc none names

%% Tongue- vs. head-centered neuronal divergence on L2 with bootstrap - Extended Data Fig. 8h

params.align_lick_num = 2;
params.bin_window = [-0.30 0.30];
params.bin_size = 0.001;
params.bin_size_back = 4;
params.bin_size_forward = 0;
params.bin_step = 1;
params.num_bins_sig = 5;%5
params.pval = 0.05;
params.bootstrap = true;
params.num_shuffles = 1000;
params.par_shuf = true;
params.ignore_first_bin = true;

[wrs_hist_tc, hist_bins_tc, ~, min_wrs_latency_sig_tc, min_wrs_latency_boot_tc] = get_contact_response_latency(sp_vars_struct_tc, selectivity_tc, params, lick_num);
[wrs_hist_hc, hist_bins_hc, ~, min_wrs_latency_sig_hc, min_wrs_latency_boot_hc] = get_contact_response_latency(sp_vars_struct_hc, selectivity_hc, params, lick_num);
[wrs_hist_both, hist_bins_both, ~, min_wrs_latency_sig_both, min_wrs_latency_boot_both] = get_contact_response_latency(sp_vars_struct_both, selectivity_both, params, lick_num);
[wrs_hist_none, hist_bins_none, ~, min_wrs_latency_sig_none, min_wrs_latency_boot_none] = get_contact_response_latency(sp_vars_struct_none, selectivity_none, params, lick_num);

med_latency_none = median(min_wrs_latency_sig_none, 'omitnan');
iqr_latency_none = prctile(min_wrs_latency_boot_none, [25 75]);
med_latency_tc = median(min_wrs_latency_sig_tc, 'omitnan');
iqr_latency_tc = prctile(min_wrs_latency_boot_tc, [25 75]);
med_latency_both = median(min_wrs_latency_sig_both, 'omitnan');
iqr_latency_both = prctile(min_wrs_latency_boot_both, [25 75]);
med_latency_hc = median(min_wrs_latency_sig_hc, 'omitnan');
iqr_latency_hc = prctile(min_wrs_latency_boot_hc, [25 75]);

% test how many of the bootstrapped medians are greater than the
% head-centric frame (p-value)
sum(min_wrs_latency_boot_none > med_latency_hc)
sum(min_wrs_latency_boot_tc > med_latency_hc)
sum(min_wrs_latency_boot_both > med_latency_hc)

figure();

min_wrs_latency_all = {min_wrs_latency_sig_none; min_wrs_latency_sig_tc; min_wrs_latency_sig_both; min_wrs_latency_sig_hc};
for i = 1:4
    min_wrs_latency_temp = min_wrs_latency_all{i};
    for j = 1:numel(min_wrs_latency_temp)
        scatter(i, min_wrs_latency_temp(j), [], [220/255 220/255 220/255], 'filled'); hold on;
    end
end

xlim([0 5]);
ylim([0 250]);
xticks([1 2 3 4])
xticklabels({'Sig, No Corr', 'Tongue-Centric', 'Tongue & Head', 'Head-Centric'})
ylabel('Divergence Latency (ms)')

scatter(1, med_latency_none, [], 'black', 'filled'); 
scatter(2, med_latency_tc, [], 'black', 'filled'); 
scatter(3, med_latency_both, [], 'black', 'filled'); 
scatter(4, med_latency_hc, [], 'black', 'filled'); 
errorbar(1, med_latency_none, med_latency_none-iqr_latency_none(1), iqr_latency_none(2) - med_latency_none, 'black');
errorbar(2, med_latency_tc, med_latency_tc-iqr_latency_tc(1), iqr_latency_tc(2) - med_latency_tc, 'black');
errorbar(3, med_latency_both, med_latency_both-iqr_latency_both(1), iqr_latency_both(2) - med_latency_both, 'black');
errorbar(4, med_latency_hc, med_latency_hc-iqr_latency_hc(1), iqr_latency_hc(2) - med_latency_hc, 'black');

%% Scatter of L2 and L4 correlation coefficient - Extended Data Fig. 8c

% store p-values and corr. coefficients of each neuron. plot if desired
corr_lick_nums = [2 4];
corr_r_L2L4 = nan(size(contact_angle_bin, 1), numel(corr_lick_nums));
corr_p_L2L4 = nan(size(contact_angle_bin, 1), numel(corr_lick_nums));
corr_coef_shuf_L2L4 = nan(size(contact_angle_bin, 1), numel(corr_lick_nums), num_shuffles);
p_shuff_L2L4 = nan(size(contact_angle_bin, 1), numel(corr_lick_nums));
r_shuff_L2L4_temp = nan(size(contact_angle_bin, 1), numel(corr_lick_nums), num_shuffles);
square_r = false;

for i = 1:size(contact_angle_bin, 1)
    for j = 1:numel(corr_lick_nums)
        [r, p] = corrcoef(-13.5:3:13.5, firing_rate_bin(i, lick_num == corr_lick_nums(j), 2:11), 'rows', 'complete');
        if square_r == true
            corr_r_L2L4(i, j) = r(1, 2).^2;
        else
            corr_r_L2L4(i, j) = r(1, 2);
        end
        corr_p_L2L4(i, j) = p(1, 2);
    
        for k = 1:num_shuffles
            temp = corrcoef(-13.5:3:13.5, firing_rate_bin_shuff(i, lick_num == corr_lick_nums(j), 2:11, k), 'rows', 'complete');
            r_shuff_L2L4_temp(i, j, k) = temp(1, 2).^2;
            corr_coef_shuf_L2L4(i, j, k) = temp(1, 2);
        end
        
        corr_95 = prctile(squeeze(r_shuff_L2L4_temp(i, j, :)), 95);
        if square_r == true
            p_shuff_L2L4(i, j) = corr_r_L2L4(i, j) > corr_95;
        else
            p_shuff_L2L4(i, j) = corr_r_L2L4(i, j).^2 > corr_95;
        end
    end
end

% by correlation values
L2_only = (p_shuff_L2L4(:, 1) == 1 & corr_p_L2L4(:, 1) < 0.05) & (p_shuff_L2L4(:, 2) == 0 | corr_p_L2L4(:, 2) >= 0.05);
L4_only = (p_shuff_L2L4(:, 2) == 1 & corr_p_L2L4(:, 2) < 0.05) & (p_shuff_L2L4(:, 1) == 0 | corr_p_L2L4(:, 1) >= 0.05);
L2_L4 = p_shuff_L2L4(:, 1) == 1 & p_shuff_L2L4(:, 2) == 1 & corr_p_L2L4(:, 1) < 0.05 & corr_p_L2L4(:, 2) < 0.05;
not_sig = (p_shuff_L2L4(:, 1) == 0 & p_shuff_L2L4(:, 2) == 0) | (corr_p_L2L4(:, 1) >= 0.05 & corr_p_L2L4(:, 2) >= 0.05);

figure('units','normalized','outerposition',[0 0 1 1]); hold on;
scatter(corr_r_L2L4(L2_only, 1), corr_r_L2L4(L2_only, 2), 'filled', 'red');
scatter(corr_r_L2L4(L4_only, 1), corr_r_L2L4(L4_only, 2), 'filled', 'black');
scatter(corr_r_L2L4(L2_L4, 1), corr_r_L2L4(L2_L4, 2), 'filled', 'blue');
scatter(corr_r_L2L4(not_sig, 1), corr_r_L2L4(not_sig, 2));

clearvars corr_lick_nums r_shuff_L2L4_temp square_r corr_95 L2_only L4_only L2_L4 not_sig

%% Plot donutchart for neurons activated by both left/right recentering - Extended Data Fig. 8f

L4_neurons = sum([selectivity{:, 2}] ~= 0 | [selectivity{:, 3}] ~= 0);
L4_both_recenter = sum([selectivity{:, 2}] == 1 & [selectivity{:, 3}] == 3);
L4_other = L4_neurons - L4_both_recenter;

figure('units','normalized','outerposition',[0 0 1 1]);
mat_version = version;
if str2double(mat_version(1:2)) >= 23
    % MATLAB 2023B AND BEYOND: Donut chart of percentage contra vs. ispi fakeout neurons
    names = categorical({strcat('Activated by both Left & Right Re-center n=', num2str(L4_both_recenter)), strcat('Not Activated by both n=', num2str(L4_other))});
    donutchart([L4_both_recenter, L4_other], names);
else
    % PRE-MATLAB 2023B: Bar chart of percentage contra vs. ispi fakeout neurons
    names = [strcat("Activated by both Left & Right Re-center n=", num2str(L4_both_recenter)), strcat("Not Activated by both n=", num2str(L4_other))];
    pie([L4_both_recenter, L4_other], names);
end

clearvars L4_neurons L4_both_recenter L4_other 

%% Plot donutchart for retuning of center L2 neurons - Extended Data Fig. 8f

center_neurons = sum([selectivity{:, 1}] == 2);
center_retuned = sum([selectivity{:, 1}] == 2 & ([selectivity{:, 2}] ~= 0 | [selectivity{:, 3}] ~= 0));
center_not_retuned = center_neurons - center_retuned;

figure('units','normalized','outerposition',[0 0 1 1]);
mat_version = version;
if str2double(mat_version(1:2)) >= 23
    % MATLAB 2023B AND BEYOND: Donut chart of percentage contra vs. ispi fakeout neurons
    names = categorical({strcat('Retuned n=', num2str(center_retuned)), strcat('Not Retuned n=', num2str(center_not_retuned))});
    donutchart([center_retuned, center_not_retuned], names);
else
    % PRE-MATLAB 2023B: Bar chart of percentage contra vs. ispi fakeout neurons
    names = [strcat("Retuned n=", num2str(center_retuned)), strcat("Not Retuned n=", num2str(center_not_retuned))];
    pie([L4_both_recenter, L4_other], names);
end

clearvars center_neurons center_retuned center_not_retuned 

%% Function List

function plot_percent_fakeout_recenter(sig_struct, sig_lick_num)

fakeout_sig = sig_struct.fakeout_sig;
recenter_sig = sig_struct.recenter_sig;

% find the number/percent significant for fakeout only. note: fakeout_sig 
% contains significance across all 3 conditions, but selectivity contains
% only left or right fakeout significant.  these numbers may differ.
num_L2_fakeout_sig = sum(fakeout_sig(:, sig_lick_num == 2) == 1);

% number and percent for recenter sig only
num_L4_recenter_sig = sum(recenter_sig(:, 1) | recenter_sig(:, 2));

% number and percent for fakeout and recenter sig
num_fakeout_recenter_sig = sum((fakeout_sig(:, sig_lick_num == 2) == 1) & (recenter_sig(:, 1) | recenter_sig(:, 2)));

% number not modulated by fakeout/recenter
num_not_fakeout_or_recenter_sig = size(fakeout_sig, 1) - (num_L2_fakeout_sig + num_L4_recenter_sig - num_fakeout_recenter_sig);

% number modulated by only fakeout
num_only_fakeout_sig = num_L2_fakeout_sig - num_fakeout_recenter_sig;

% number modulated by only recenter
num_only_recenter_sig = num_L4_recenter_sig - num_fakeout_recenter_sig;

% combine data for each condition
data = [num_only_fakeout_sig num_only_recenter_sig num_fakeout_recenter_sig num_not_fakeout_or_recenter_sig];

% create labels for each condition
fakeout_only_label = strcat("Fakeout Only (n=", num2str(num_only_fakeout_sig), ")");
recenter_only_label = strcat("Recenter Only (n=", num2str(num_only_recenter_sig), ")");
fakeout_recenter_label = strcat("Fakeout & Recenter (n=", num2str(num_fakeout_recenter_sig), ")");
not_modulated_label = strcat("Not Modulated (n=", num2str(num_not_fakeout_or_recenter_sig), ")");
labels = [fakeout_only_label, recenter_only_label, fakeout_recenter_label, not_modulated_label];

mat_version = version;
if str2double(mat_version(1:2)) >= 23
    % MATLAB 2023B AND BEYOND:  Donut chart of percent modulated by fakeout/recenter only, fakeout + recenter or not modulated
    figure('units','normalized','outerposition',[0 0 1 1]);
    donutchart(data, labels);
else
    % PRE-MATLAB 2023B:  Pie chart of percent modulated by fakeout/recenter only, fakeout + recenter or not modulated
    figure(); pie(data, labels);
end
end

function [fiducial_session_date, fiducial_session_mouse] = parse_SC_Ephys_mouse_date_fiducial_list(fiducial_list)

fiducial_session_date = cell(size(fiducial_list, 1), 1);
fiducial_session_mouse = cell(size(fiducial_list, 1), 1);
for i = 1:size(fiducial_list, 1)
    fiducial_session_date_temp = regexp(fiducial_list{i, 1}, '(?<=Masks\\).*', 'match');
    fiducial_session_date_temp = split(fiducial_session_date_temp, "_");
    fiducial_session_date{i} = fiducial_session_date_temp{1};

    fiducial_session_mouse_temp = regexp(fiducial_list{i, 1}, '(?<=SC_Ephys\\).*', 'match');
    if isempty(fiducial_session_mouse_temp)
        fiducial_session_mouse_temp = regexp(fiducial_list{i, 1}, '(?<=Anesthesia\\).*', 'match');
        if isempty(fiducial_session_mouse_temp)
            fiducial_session_mouse_temp = regexp(fiducial_list{i, 1}, '(?<=SpXII_Ephys\\).*', 'match');
        end
    end
    fiducial_session_mouse_temp = split(fiducial_session_mouse_temp, "\");
    fiducial_session_mouse{i} = fiducial_session_mouse_temp{1};
end
end

function plot_same_nick_site_fakeout_recenter(sig_struct)

selectivity = sig_struct.selectivity_cent;

% number of neurons who were activated for the same type of fakeout nick
% and recenter nick
num_same_left = sum(([selectivity{:, 1}] == 1 | [selectivity{:, 1}] == 4) & [selectivity{:, 2}] == 1);
num_same_right = sum(([selectivity{:, 1}] == 3 | [selectivity{:, 1}] == 5) & [selectivity{:, 3}] == 3);

num_same_fakeout_recenter_sig = num_same_left + num_same_right;
total_fakeout_and_recenter_modulated = sum(([selectivity{:, 1}] > 0) & ([selectivity{:, 2}] ~= 0 | [selectivity{:, 3}] ~= 0));
number_other_sig = (total_fakeout_and_recenter_modulated - num_same_fakeout_recenter_sig);

same_label = strcat("Activated by Same Nick on Fakeout and Recenter (n = ", num2str(num_same_fakeout_recenter_sig), ")");
other_label = strcat("Other (n = ", num2str(number_other_sig), ")");
labels_activated_fakeout_recenter = categorical([same_label, other_label]);

mat_version = version;
if str2double(mat_version(1:2)) >= 23
    % MATLAB 2023B AND BEYOND: Donut chart of percentage activated by the same nick side
    figure('units','normalized','outerposition',[0 0 1 1]);
    donutchart([num_same_fakeout_recenter_sig, number_other_sig], labels_activated_fakeout_recenter);
else
    % PRE-MATLAB 2023B: Bar chart of percentage contra vs. ispi fakeout neurons
    figure(); bar(labels_activated_fakeout_recenter, [num_same_fakeout_recenter_sig, number_other_sig], 0.5);
    ylabel('Number of Significant Fakeout Neurons');
    ylim([0 100]);
end
end

function [sp_session_date, sp_session_mouse] = parse_SC_Ephys_mouse_date_selectivity(selectivity, no_recenter)

sp_session_date = cell(size(selectivity, 1), 1);
sp_session_mouse = cell(size(selectivity, 1), 1);
for i = 1:size(selectivity, 1)
    if no_recenter == false
        sp_session_date_temp = split(selectivity{i, 5}, "_");
       
        sp_session_date{i} = sp_session_date_temp{1};
        
        sp_session_mouse_start_ind = regexp(selectivity{i, 5}, strcat('(?<=', sp_session_date_temp{1}, '_).*'));
        [~, sp_session_mouse_end_ind] = regexp(selectivity{i, 5}, '\w*(?=_recentering)');
        sp_session_mouse{i} = selectivity{i, 5}(sp_session_mouse_start_ind:sp_session_mouse_end_ind); 

        % I fucked up  named 230806 230805 for SC_Ephys_4...temp fix here.
        if strcmp(sp_session_date_temp{1}, '230805') & strcmp(sp_session_mouse{i} , 'SC_Ephys_4')
            sp_session_date{i} = '230806';
        end
        
    else 
        sp_session_date_temp = split(selectivity{i, 3}, "_");
       
        sp_session_date{i} = sp_session_date_temp{1};
        
        sp_session_mouse_start_ind = regexp(selectivity{i, 3}, strcat('(?<=', sp_session_date_temp{1}, '_).*'));
        [~, sp_session_mouse_end_ind] = regexp(selectivity{i, 3}, '\w*(?=_recentering)');
        sp_session_mouse{i} = selectivity{i, 3}(sp_session_mouse_start_ind:sp_session_mouse_end_ind);  

        % I fucked up  named 230806 230805 for SC_Ephys_4...temp fix here.
        if strcmp(sp_session_date_temp{1}, '230805') & strcmp(sp_session_mouse{i} , 'SC_Ephys_4')
            sp_session_date{i} = '230806';
        end

        % we didn't name Anesthesia mice properly, this is the workaround
        if contains(selectivity{i, 3},'230131_ALM') || contains(selectivity{i, 3},'230201_ALM') || contains(selectivity{i, 3},'230202_ALM') || contains(selectivity{i, 3},'230204_ALM') || contains(selectivity{i, 3},'230205_ALM') || contains(selectivity{i, 3},'230206_ALM')
            sp_session_mouse{i} = 'Anesthesia_2';
        elseif contains(selectivity{i, 3},'ALM_fakeout2D') || contains(selectivity{i, 3},'SC_fakeout2D') 
            sp_session_mouse{i} = 'Anesthesia_1';
        elseif contains(selectivity{i, 3},'fakeout2D_ephys_SC')  
            sp_session_mouse{i} = 'SpXII_Ephys_2';
        end
    end
    
end
end

function [behavior_bin, firing_rate_bin, firing_rate_bin_shuff, num_trials_bin, mean_fr_shuff, mean_fr] = bin_behavior_and_firing_rate(behavior_sig, firing_rate_sig, bin, min_trial_num, num_shuffles)

if ~exist('num_shuffles', 'var')
    shuffle_enable = 0;
else
    shuffle_enable = 1;
end

behavior_bin = [];
firing_rate_bin = [];
behavior_bin_shuff = [];
firing_rate_bin_shuff = [];
num_trials_bin = [];
mean_fr = [];
mean_fr_shuff = [];
for i = 1:size(behavior_sig, 2)
    for j = 1:size(behavior_sig{i}, 1)
        n = 1;
        firing_rate_all = {};
        for k = 1:(numel(bin) - 1)
            lat_disp_change_temp = behavior_sig{i}(j, behavior_sig{i}(j, :) >= bin(k) & behavior_sig{i}(j, :) < bin(k + 1));
            if numel(lat_disp_change_temp) >= min_trial_num
                num_trials_bin(i, j, k) = numel(lat_disp_change_temp);
                behavior_bin(i, j, k) = mean(lat_disp_change_temp);
                
                firing_rate_sort = firing_rate_sig{i}(j, behavior_sig{i}(j, :) >= bin(k) & behavior_sig{i}(j, :) < bin(k + 1));
                firing_rate_bin(i, j, k) = mean(firing_rate_sort);
                firing_rate_all{n} = firing_rate_sort;
                n = n + 1;

            else
                behavior_bin(i, j, k) = NaN;
                firing_rate_bin(i, j, k) = NaN;
                num_trials_bin(i, j, k) = NaN;
            end
        end

        mean_fr(i, j) = mean([firing_rate_all{:}]);

        for m = 1:num_shuffles
            firing_rate_ind = randperm(size(firing_rate_sig{i}, 2));
            firing_rate_sig_shuff = firing_rate_sig{i}(j, firing_rate_ind);
            n = 1;
            firing_rate_all_shuff = {};
            for k = 1:(numel(bin) - 1)
                lat_disp_change_temp = behavior_sig{i}(j, behavior_sig{i}(j, :) >= bin(k) & behavior_sig{i}(j, :) < bin(k + 1));
                if numel(lat_disp_change_temp) >= min_trial_num
                    firing_rate_sort_shuff = firing_rate_sig_shuff(1, behavior_sig{i}(j, :) >= bin(k) & behavior_sig{i}(j, :) < bin(k + 1));
                    firing_rate_bin_shuff(i, j, k, m) = mean(firing_rate_sort_shuff);
                    firing_rate_all_shuff{n} = firing_rate_sort_shuff;
                    n = n + 1;

                else
                    firing_rate_bin_shuff(i, j, k, :) = NaN;
                end
            end
            mean_fr_shuff(i, j, m) = mean([firing_rate_all_shuff{:}]);
        end
    end
end
end

function site_sel_sig = plot_site_selectivity(selectivity_bin, selectivity_rand_bin, num_shuffles)

center_bin = selectivity_bin{1};
contra_bin = selectivity_bin{2};
ipsi_bin = selectivity_bin{3};
center_contra_bin = selectivity_bin{4};
center_ipsi_bin = selectivity_bin{5};

center_bin_rand = selectivity_rand_bin{1};
contra_bin_rand = selectivity_rand_bin{2};
ipsi_bin_rand = selectivity_rand_bin{3};
center_contra_bin_rand = selectivity_rand_bin{4};
center_ipsi_bin_rand = selectivity_rand_bin{5};

center_iqr = prctile(center_bin_rand', [25 75]);
contra_iqr = prctile(contra_bin_rand', [25 75]);
ipsi_iqr = prctile(ipsi_bin_rand', [25 75]);
center_contra_iqr = prctile(center_contra_bin_rand', [25 75]);
center_ipsi_iqr = prctile(center_ipsi_bin_rand', [25 75]);

for i = 1:numel(center_bin)
    y_temp(i, 1:5) = [contra_bin(i) center_contra_bin(i) center_bin(i) center_ipsi_bin(i) ipsi_bin(i)];
    neg(i, 1:5) = y_temp(i, 1:5) - [contra_iqr(1, i) center_contra_iqr(1, i) center_iqr(1, i) center_ipsi_iqr(1, i) ipsi_iqr(1, i)];
    pos(i, 1:5) = [contra_iqr(2, i) center_contra_iqr(2, i) center_iqr(2, i) center_ipsi_iqr(2, i) ipsi_iqr(2, i)] - y_temp(i, 1:5);
end

figure('units','normalized','outerposition',[0 0 1 1]);

b = bar(y_temp); hold on;
[ngroups, nbars] = size(y_temp);

x = nan(nbars, ngroups);
for i = 1:nbars
    x(i, :) = b(i).XEndPoints;
end

errorbar(x', y_temp, neg, pos, 'k', 'linestyle', 'none')

ylim([0. 0.6]);
yticks([0 0.6]);
xticks(1:5);
xticklabels({'-3.0', '-3.2', '-3.4', '-3.6', '-3.8'});
xlabel('latSC AP Recording Site (mm)')
ylabel('Probability')
legend('Contra', 'Center+Contra', 'Center', 'Center+Ipsi', 'Ipsi')

contra_center_contra = sum(contra_bin_rand < center_contra_bin_rand, 2)/num_shuffles;
contra_center = sum(contra_bin_rand < center_bin_rand, 2)/num_shuffles;
contra_center_ipsi = sum(contra_bin_rand < center_ipsi_bin_rand, 2)/num_shuffles;
contra_ipsi = sum(contra_bin_rand < ipsi_bin_rand, 2)/num_shuffles;

center_contra_center = sum(center_contra_bin_rand < center_bin_rand, 2)/num_shuffles;
center_contra_center_ipsi = sum(center_contra_bin_rand < center_ipsi_bin_rand, 2)/num_shuffles;
center_contra_ipsi = sum(center_contra_bin_rand < ipsi_bin_rand, 2)/num_shuffles;

center_center_ipsi = sum(center_bin_rand < center_ipsi_bin_rand, 2)/num_shuffles;
center_ipsi = sum(center_bin_rand < ipsi_bin_rand, 2)/num_shuffles;

center_ipsi_ipsi = sum(center_ipsi_bin_rand < ipsi_bin_rand, 2)/num_shuffles;

site_sel_sig = table(contra_center_contra, contra_center, contra_center_ipsi, contra_ipsi, center_contra_center, center_contra_center_ipsi, center_contra_ipsi, center_center_ipsi, center_ipsi, center_ipsi_ipsi);

end
