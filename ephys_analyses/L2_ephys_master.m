%% Load full latSC dataset. All neurons are used for L2 analyses. 

load('latSC_ephys.mat')

%% Remove neurons with trials < min_trial_num & filter behavior vars

clearvars -except sp_struct

% this boolean indicates whether we are not (true) or are (false) analyzing
% re-centering sessions. In functions below, this is used to perform
% specific analyses/filtering. 
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

%% Test for significance in fakeout and recenter conditions

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

%% Plot percentage significant for fakeout

plot_percent_fakeout(sig_struct)

%% Plot percentage contra vs. ipsi selective - Fig. 3c

% note that, center selectivity here is defined as no modulation between
% left vs. right, but modulation with center only. Neurons modulated most
% by center, but modulated also between left and right are defined as left
% or right selective, not center. 

plot_contra_ipsi_fakeout(sig_struct)

%% Plot selectivity heatmaps

params.align_lick_num = 2;
params.smooth_bin_window = [-0.35 0.35];
params.real_window = [-0.3 0.3];
params.bin_size = 0.01;
params.bin_size_back = 3;
params.bin_size_forward = 0;

create_selectivity_heatmaps(sp_vars_struct, sig_struct, params, lick_num)

%% Calculate contact 'error' (e.g., contact location)

% change with the appropriate file path to Accessory Files
load('Z:\bsi8\Ito, Gao et. al (2024)\Code Review\Accessory Files\L2_ephys_fiducial_list.mat')

% calculate other kinematic variables
contact_vars_struct = calculate_contact_error_vars(sp_vars_struct, fiducial_xy);

clearvars fiducial_xy

%% Single Neuron Latency to Divergence for L2 - Fig. 3g

params.align_lick_num = 2;
params.bin_window = [-0.30 0.30];
params.bin_size = 0.001;
params.bin_size_back = 4;
params.bin_size_forward = 0;
params.bin_step = 1;
params.num_bins_sig = 3;
params.pval = 0.05/(3*params.num_bins_sig);
params.bootstrap = false;
params.num_shuffles = 1000;
params.par_shuf = true;
params.ignore_first_bin = true;

[wrs_hist, hist_bins, ~, min_wrs_latency_sig, ~] = get_contact_response_latency(sp_vars_struct, selectivity, params, lick_num);

median_wrs_latency = median(min_wrs_latency_sig, "omitnan");
mode_wrs_latency = mode(min_wrs_latency_sig);

% plot the distribution and the median
figure('units','normalized','outerposition',[0 0 1 1]);
stairs(hist_bins(1:numel(hist_bins) - 1) - 1, wrs_hist./(sum(wrs_hist)));
hold on;
plot([median_wrs_latency median_wrs_latency], [0 0.15])
ylim([0 0.15])
xlim([0 250])
xlabel('Time from L2 Onset (ms');
xlabel('Time from L2 Onset (ms)');
ylabel('Percentage of Neurons with L2 Selectivity');
yticks([0 0.15])
xticks([0 50 100 150 200 250])

clearvars params

%% PCA Divergence - Fig. 3f

num_shuffles = 1000;
pca_fakeout(sp_vars_struct, lick_num, num_shuffles)
view(-135, 45)

%% Calculate binned contact angle and firing rate correlations

contact_ang_sig = {contact_vars_struct(:).contact_ang};
firing_rate_sig = sig_struct.firing_rate(:, :, :);

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

clearvars contact_ang_sig_hc contact_ang_sig min_trial_num ang_min ang_max ang_bin_width ang_bin x

%% Calculate Skaggs Spatial Information (bits/spk) - Fig. 3d

% calculate spatial info
corr_lick_num = 2;
[contact_info_all, contact_info_shuff, contact_info_bin, contact_info_null_all] = skaggs_spatial_info(firing_rate_sig, firing_rate_bin, firing_rate_bin_shuff, num_trials_bin, mean_fr_shuff, mean_fr, corr_lick_num);

% calculate number of neurons with significant spatial information, which
% we define as greater than the 95th percentile of the null/shufffled dist.
sig_info = [];
for i = 1:numel(contact_info_all)
    prctile_temp = prctile(contact_info_shuff(i, :), 95);
    sig_info(i) = contact_info_all(i) > prctile_temp & contact_info_all(i) > contact_info_null_all(i);
end

% filter neurons with mean firing rate less than the minimum. We chose to
% make this a higher threshold than in significance testing to filter out
% very low firing rate neurons that could give rise to spurious
% correlations and contact/spatial info. 
mean_fr_min = 3;
mean_fr_ind = sum(mean_fr > mean_fr_min, 2) > 0;
mean_info = median(contact_info_all(sig_info == 1 & mean_fr_ind'));
mean_info_not_sig = median(contact_info_all(sig_info == 0 & mean_fr_ind'));

% plot the results
total_sig_info = sum(sig_info == 1 & mean_fr_ind');
total_not_sig = sum(sig_info == 0 & mean_fr_ind');

figure('units','normalized','outerposition',[0 0 1 1]);
mat_version = version;
if str2double(mat_version(1:2)) >= 23
    % MATLAB 2023B AND BEYOND: Donut chart of percentage contra vs. ispi fakeout neurons
    labels_temp = categorical({strcat('Sig Info n=', num2str(total_sig_info)), strcat('Not Sig Info n=', num2str(total_not_sig))});
    donutchart([total_sig_info, total_not_sig], labels_temp);
else
    % PRE-MATLAB 2023B: Bar chart of percentage contra vs. ispi fakeout neurons
    labels_temp = [strcat("Sig Info n=", num2str(total_sig_info)), strcat("Not Sig Info n=", num2str(total_not_sig))];
    pie([total_sig_info, total_not_sig], labels_temp);
end

figure('units','normalized','outerposition',[0 0 1 1]);
scatter(2, mean_info); hold on;
scatter(4, mean_info_not_sig);
iqr_spatial_info = prctile(contact_info_all(sig_info == 1 & mean_fr_ind'), [25 75]);
errorbar(2, mean_info, mean_info-iqr_spatial_info(1), iqr_spatial_info(2) - mean_info);
iqr_spatial_info = prctile(contact_info_all(sig_info == 0 & mean_fr_ind'), [25 75]);
errorbar(4, mean_info_not_sig, mean_info_not_sig-iqr_spatial_info(1), iqr_spatial_info(2) - mean_info_not_sig);
xlim([1 5]); 
ylim([0 0.4]);
xticks(2);
yticks([0 0.4]);
ylabel('Contact Location Information (bits/spike)')

clearvars prctile_temp mean_fr_min mean_info mean_info_not_sig total_sig_info total_not_sig iqr_spatial_info i

%% Calculate contact angle and firing rate correlation. Plot single neurons for Fig. 3a and Extended Data Fig. 6a.

% store p-values and corr. coefficients of each neuron. plot if desired
corr_lick_num = 2;
num_shuffles = 1000;
corr_r = nan(1, size(contact_angle_bin, 1));
corr_p = nan(1, size(contact_angle_bin, 1));
corr_coef_shuf = nan(1, size(contact_angle_bin, 1), num_shuffles);
p_shuff = nan(1, size(contact_angle_bin, 1));
r2_shuff_temp = nan(1, size(contact_angle_bin, 1), num_shuffles);

for i = 1:size(contact_angle_bin, 1)
    [r, p] = corrcoef(-13.5:3:13.5, firing_rate_bin(i, lick_num == corr_lick_num, 2:11), 'rows', 'complete');
    corr_r(1, i) = r(1, 2).^2;
    corr_p(1, i) = p(1, 2);

    for j = 1:num_shuffles
        temp = corrcoef(-13.5:3:13.5, firing_rate_bin_shuff(i, lick_num == corr_lick_num, 2:11, j), 'rows', 'complete');
        r2_shuff_temp(1, i, j) = temp(1, 2).^2;
        corr_coef_shuf(1, i, j) = temp(1, 2);
    end
    
    corr_95 = prctile(r2_shuff_temp(1, i, :), 95);
    p_shuff(1, i) = corr_r(1, i) > corr_95;

end

% % plot a single neuron's correlation
% neuron = 1;
% figure('units','normalized','outerposition',[0 0 1 1]); 
% scatter(-13.5:3:13.5, squeeze(firing_rate_bin(neuron, lick_num == corr_lick_num, 2:11)));
% xlim([-15 15])

clearvars plot_neuron square_r x contact_angle_bin_temp firing_rate_bin_temp firing_rate_bin_shuff_temp i...
    r p j temp corr_95 corr_lick_num

%% Plotting correlations - Fig. 3e

corr_r_temp = corr_r(:, mean_fr_ind);
corr_p_temp = corr_p(:, mean_fr_ind);
p_shuff_temp = p_shuff(:, mean_fr_ind);

% Calculate percentage significant correlation relative to shuffles for all
sig_corr_shuff = sum(p_shuff_temp(1, :) == 1 & corr_p_temp(1, :) < 0.05)/size(p_shuff_temp, 2);
not_sig_corr_shuff = sum(p_shuff_temp(1, :) ~= 1 | corr_p_temp(1, :) >= 0.05)/size(p_shuff_temp, 2);

figure();
mat_version = version;
if str2double(mat_version(1:2)) >= 23
    % MATLAB 2023B AND BEYOND: Donut chart of percentage contra vs. ispi fakeout neurons
    labels_temp = categorical({'Sig Corr', 'Not Sig Corr'});
    donutchart([sig_corr_shuff, not_sig_corr_shuff], labels_temp);
else
    % PRE-MATLAB 2023B: Bar chart of percentage contra vs. ispi fakeout neurons
    labels_temp = ["Sig Corr", "Not Sig Info"];
    pie([sig_corr_shuff, not_sig_corr_shuff], labels_temp);
end

% Plot median + IQR of R2 for significantly correlated neurons
sig_corr_r = corr_r_temp(1, p_shuff_temp(1, :) == 1 | corr_p_temp(1, :) < 0.05);
median_corr = nanmedian(sig_corr_r);
not_sig_corr_r = corr_r_temp(1, p_shuff_temp(1, :) == 0 | corr_p_temp(1, :) > 0.05);
median_corr_not_sig = nanmedian(not_sig_corr_r);
iqr_corr = prctile(sig_corr_r, [25 75]);
iqr_not_corr = prctile(not_sig_corr_r, [25 75]);
figure('units','normalized','outerposition',[0 0 1 1]);
scatter(2, median_corr); hold on;
scatter(4, median_corr_not_sig)
errorbar(2, median_corr, median_corr - iqr_corr(1), iqr_corr(2) - median_corr);
errorbar(4, median_corr_not_sig, median_corr_not_sig - iqr_not_corr(1), iqr_not_corr(2) - median_corr_not_sig);

xlim([1 5])
ylim([0 1])

clearvars corr_r_temp corr_p_temp p_shuff_temp corr_alpha median_corr not_sig_corr_r median_corr_not_sig iqr_corr iqr_not_corr

%% Cross-correlation of lick cycle with firing rates - Extended Data Fig. 9b - d

clearvars -except sp_struct

% this boolean indicates whether we are not (true) or are (false) analyzing
% re-centering sessions. In functions below, this is used to perform
% specific analyses/filtering. 
no_recenter = true;

% define lick numbers to be analyzed and params used for filtering trials.
% These analyses are focused on L2 and L3. We require L1 to L3 to exist
% (lick_exist_lick_num), to make spout contact on each lick
% (sp_contact_lick_num) only once (sp_contact2_absent_lick_num), and have
% contact tracking data (contact_centroid_lick_num).  
lick_num = [1 2 3 4 5];
params.lick_exist_lick_num = [1 2 3 4 5];
params.sp_contact_lick_num = [1 2 3 4 ];
params.sp_contact2_absent_lick_num = [1 2 3 4 5];
params.contact_centroid_lick_num = [1 2 3 4 5];

% exclude neurons with trials < min_trial_num & firing rate < min_fr_thresh,
% for trials that pass the behavioral params set above. 
min_trial_num = 10;
min_fr_thresh = 0.1;
[sp_vars_struct, sp_struct_filt] = filter_sp_neurons_and_vars_recenter(sp_struct, min_trial_num, min_fr_thresh, lick_num, params, no_recenter);

clearvars min_trial_num min_fr_thresh params

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

% This code analyzes cross-correlation between lick cycles (quantified with
% tongue volume) and neurons' firing rates. This section requires having 
% sp_vars_struct + sig_struct, as derived in sig_test_firing_rates_recenter_2D_script

% ***Note: Significant cross-correlation is determined from bootstrapping 
% and shuffling. Therefore, each run of this code may yield slightly
% different significance outcome for a very small percentage of neurons.
% However, none of the differences should affect the conclusions drawn from
% this analysis

neuron_ID = []; % leave empty if all neurons need to be analyzed
start_lick = 4; % first lick to be analyzed/the lick whose protrusion onset everyting  is aligned to
analysis_window = [50 300]; % when to analyze, relative to start lick onset, in units if ms
lag_window = [-100 100]; % the lag window to use in cross-correlation
trial_types = [1,2,3]; % if all ipsi&center&contra trials need to be analyzed
bin_size = 10;
num_shuffle = 1000; % for shuffles and bootstrapps
plot_individual = 0; % if individual plots of each neuron are needed
[selectivity_aiming,r_aiming,lag_aiming] = firing_rate_lick_cycle_xcorr(sp_vars_struct,sig_struct,neuron_ID,start_lick,analysis_window,lag_window,trial_types,bin_size,num_shuffle,plot_individual);

plot_xcorr_summary(selectivity_aiming,r_aiming,lag_aiming)

% function below collects data needed for a three-circle venn diagram,showing neurons with significant xcorr for each trial type
venn_data = get_venn_data(selectivity_aiming); % legened: ipsi, center, contra,ipsi&center, ipsi&contra, center&contra, all three, none

%% Function List

function plot_percent_fakeout(sig_struct)

fakeout_sig = [sig_struct.selectivity_cent{:, 1}];

% find the number/percent significant for fakeout only. note: fakeout_sig 
% contains significance across all 3 conditions, but selectivity contains
% only left or right fakeout significant.  these numbers may differ.
num_L2_fakeout_sig = sum(fakeout_sig ~= 0);

% number not modulated by fakeout/recenter
num_not_fakeout_sig = numel(fakeout_sig) - (num_L2_fakeout_sig);

% combine data for each condition
data = [num_L2_fakeout_sig num_not_fakeout_sig];

% create labels for each condition
fakeout_only_label = strcat("Fakeout (n=", num2str(num_L2_fakeout_sig), ")");
not_modulated_label = strcat("Not Modulated (n=", num2str(num_not_fakeout_sig), ")");
labels = [fakeout_only_label, not_modulated_label];

figure('units','normalized','outerposition',[0 0 1 1]);
mat_version = version;
if str2double(mat_version(1:2)) >= 23
    % MATLAB 2023B AND BEYOND:  Donut chart of percent modulated by fakeout/recenter only, fakeout + recenter or not modulated
    donutchart(data, labels);
else
    % PRE-MATLAB 2023B:  Pie chart of percent modulated by fakeout/recenter only, fakeout + recenter or not modulated
    pie(data, labels);
end
end

function plot_contra_ipsi_fakeout(sig_struct)

selectivity = sig_struct.selectivity_cent;

% calculate number of neurons L2 fakeout significant
num_L2_fakeout_sig = sum([selectivity{:, 1}] ~= 0);

% get whether the significant left neurons were recorded from right hemi
left_session_name = selectivity([selectivity{:, 1}] == 1, 3);
sig_left_contra = cellfun(@(x) contains(x, 'right'), left_session_name);

center_left_session_name = selectivity([selectivity{:, 1}] == 4, 3);
sig_center_left_contra = cellfun(@(x) contains(x, 'right'), center_left_session_name);

% get whether the significant left neurons were recorded from right hemi
right_session_name = selectivity([selectivity{:, 1}] == 3, 3);
sig_right_contra = cellfun(@(x) contains(x, 'left'), right_session_name);

center_right_session_name = selectivity([selectivity{:, 1}] == 5, 3);
sig_center_right_contra = cellfun(@(x) contains(x, 'left'), center_right_session_name);

sig_contra = [sig_left_contra; sig_right_contra];
sig_center_contra = [sig_center_left_contra; sig_center_right_contra];
num_sig_contra = sum(sig_contra);
num_sig_ipsi = numel(sig_contra) - num_sig_contra;
num_sig_center = sum([selectivity{:, 1}] == 2);
num_sig_center_contra = sum(sig_center_contra);
num_sig_center_ipsi = numel(sig_center_contra) - sum(sig_center_contra);

contra_label = strcat("Contra Selective (n = ", num2str(num_sig_contra), ")");
ipsi_label = strcat("Ipsi Selective (n = ", num2str(num_sig_ipsi), ")");
center_label = strcat("Center Selective (n = ", num2str(num_sig_center), ")");
center_contra_label = strcat("Center Contra Selective (n = ", num2str(num_sig_center_contra), ")");
center_ipsi_label = strcat("Center Ipsi Selective (n = ", num2str(num_sig_center_ipsi), ")");

labels_fakeout = categorical([contra_label, center_contra_label, center_label, ipsi_label, center_ipsi_label]);
number_contra_ipsi_fakeout = [num_sig_contra, num_sig_center_contra, num_sig_center, num_sig_ipsi, num_sig_center_ipsi];

figure('units','normalized','outerposition',[0 0 1 1]);
mat_version = version;
if str2double(mat_version(1:2)) >= 23
    % MATLAB 2023B AND BEYOND: Donut chart of percentage contra vs. ispi fakeout neurons
    donutchart(number_contra_ipsi_fakeout, labels_fakeout);
else
    % PRE-MATLAB 2023B: Bar chart of percentage contra vs. ispi fakeout neurons
    bar(labels_fakeout, number_contra_ipsi_fakeout, 0.5);
    ylabel('Number of Significant Fakeout Neurons');
    ylim([0 100]);
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

function [selectivity_aiming,r_aiming,lag_aiming]= firing_rate_lick_cycle_xcorr(sp_vars_struct,sig_struct,neuron_ID,start_lick,analysis_window,lag_window,trial_types,bin_size,num_shuffle,plot_individual)
% initialize variables based on preset parameters
if isempty(neuron_ID)
    neuron_ID = [1:size(sp_vars_struct,2)];
end
selectivity_aiming = zeros(max(neuron_ID),3);
lag_aiming = nan(max(neuron_ID),3);
r_aiming = nan(max(neuron_ID),3);
r_shuff_aiming = nan(max(neuron_ID),3);
session_name = sig_struct.selectivity(:,end);
lag_range = lag_window(2)-lag_window(1)+1;
lag_window_binned = [lag_window(1)/bin_size:1:lag_window(2)/bin_size];
plot_fr = nan(max(neuron_ID),max(trial_types),1300/bin_size);
plot_lr = nan(max(neuron_ID),max(trial_types),1300/bin_size);
plot_xcorr = nan(max(neuron_ID),max(trial_types),1+((lag_range-1)/bin_size));

for n = neuron_ID
    sp_vars_temp = sp_vars_struct(n);
    selectivity_temp = zeros(1,3);

    for i = trial_types % separate by trial type
        % no recentering trials are used
        trial_ind = [sp_vars_temp.fakeout_trial_filt == i & sp_vars_temp.recenter_trial_filt == 0];
        trial_ind = find(trial_ind == 1);
        fr_all = nan(numel(sp_vars_temp.sp_times_filt),1300);
        lr_all = nan(numel(sp_vars_temp.sp_times_filt),1300);
                
        for j = trial_ind
            % collect lick rate (volume) and firing rate till the end of the analysis window
            train_offset = sp_vars_temp.prot_filt(start_lick,j) + analysis_window(2);
            volume_temp = sp_vars_temp.volume_filt(:,j);
            if train_offset >1300
                train_offset = 1300;
            end
            lr_train = zeros(1,train_offset);
            for k = 2:numel(volume_temp)
                if sp_vars_temp.prot_filt(k,j) > 0
                    lr_train(sp_vars_temp.prot_filt(k,j):sp_vars_temp.prot_filt(k,j)+numel(volume_temp{k})-1) = volume_temp{k};
                end
            end            
            sp_times_temp = sp_vars_temp.sp_times_filt{j};
            fr_train = zeros(1,train_offset);
            for k = 1:numel(sp_times_temp) 
                if ceil(sp_times_temp(k)) > 0 && ceil(sp_times_temp(k)) <= numel(fr_train)
                    fr_train(ceil(sp_times_temp(k))) = 1;
                end
            end
            % trim lick rate and firing rate at the start of the analysis window           
            train_onset = sp_vars_temp.prot_filt(start_lick,j) - analysis_window(1);
            if train_onset <= 0
                fr_all(j,:) = nan(1,numel(fr_all(j,:)));
                lr_all(j,:) = nan(1,numel(fr_all(j,:)));
            else
                fr_train_temp = fr_train(train_onset:end);
                fr_all(j,1:numel(fr_train_temp)) = fr_train_temp;
                lr_train_temp = lr_train(train_onset:end);
                lr_all(j,1:numel(lr_train_temp)) = lr_train_temp;
            end 
        end
        
        if mean(mean(fr_all,2,'omitnan'),'omitnan')*1000 < 5 % trials with firing rate lower than 5hz are excluded
            continue
        end

        % bin firing rate and lick rate based on the preset bin size
        bin_num = floor(size(fr_all,2)/bin_size);
        fr_train_bin = nan(size(fr_all,1),bin_num);
        lr_train_bin = nan(size(fr_all,1),bin_num);        
        for j = 1:size(fr_all,1)
            fr_train_bin_temp = nan(1,bin_num);
            lr_train_bin_temp = nan(1,bin_num);
            for k = 1:bin_num
                fr_train_bin_temp(k) = mean(fr_all(j,(k-1)*bin_size + 1 : k*bin_size));
                lr_train_bin_temp(k) = mean(lr_all(j,(k-1)*bin_size + 1 : k*bin_size));
            end
            fr_train_bin(j,:) = fr_train_bin_temp;
            lr_train_bin(j,:) = lr_train_bin_temp;
        end
        
        % average firing rate and lick rate across trials
        fr_train_xcross = nan(1,size(fr_train_bin,2));
        lr_train_xcross = nan(1,size(fr_train_bin,2));
        for j = 1:size(fr_train_bin,2)
            if sum(~isnan(fr_train_bin(:,j))) >= 10 % there should be a minimum of 10 trials
                fr_train_xcross(1,j) = mean(fr_train_bin(:,j),'omitnan')*1000;
                lr_train_xcross(1,j) = mean(lr_train_bin(:,j),'omitnan');
            end
        end
        % perform mean substraction to prepare for xcorr
        xcorr_fr = fr_train_xcross-mean(fr_train_xcross,'omitnan');
        xcorr_fr = xcorr_fr(~isnan(xcorr_fr));
        xcorr_lr = lr_train_xcross-mean(lr_train_xcross,'omitnan');
        xcorr_lr = xcorr_lr(~isnan(xcorr_lr));
        if plot_individual == 1
            plot_fr(n,i,1:numel(fr_train_xcross)) = fr_train_xcross-mean(fr_train_xcross,'omitnan');
            plot_lr(n,i,1:numel(fr_train_xcross)) = lr_train_xcross-mean(lr_train_xcross,'omitnan');
        end
        
        % run xcorr, get the maximum r value over the lag window
        [r,lags] = xcorr(xcorr_fr,xcorr_lr,'normalized');
        r = r([lags>=lag_window(1)/bin_size & lags<=lag_window(2)/bin_size]);
        r_actual = unique(max(r));
        if isempty(r_actual)
            continue
        end        
        plot_xcorr(n,i,:) = r;

        % create a distribution of max_r by randomly shuffling firing rate within each trial
        r_shuff = nan(1,num_shuffle);
        r_shuff_all = nan(num_shuffle,numel(r));
        fr_train_bin_shuffle = fr_train_bin(trial_ind,:);
        lr_train_bin_shuffle = lr_train_bin(trial_ind,:);
        for j = 1:num_shuffle
            fr_temp = nan(size(fr_train_bin_shuffle,1),size(fr_train_bin_shuffle,2));
            for k = 1:size(fr_train_bin_shuffle,1)
                rand_ind = randperm(find(~isnan(fr_train_bin_shuffle(k,:)),1,'last'));
                fr_temp(k,1:max(rand_ind)) = fr_train_bin_shuffle(k,rand_ind);
            end
            [shuff_fr,shuff_lr] = get_trains_for_xcross(fr_temp, lr_train_bin_shuffle);
            [r_temp,lags] = xcorr(shuff_fr,shuff_lr,'normalized');
            r_shuff_all(j,:) = r_temp([lags>=lag_window(1)/bin_size & lags<=lag_window(2)/bin_size]);
            r_shuff(j) = unique(max(r_shuff_all(j,:)));
        end

        % create a distribution of natural variability of max_r by bootstrapping
        r_BS = nan(1,num_shuffle);
        r_BS_all = nan(num_shuffle,numel(r));
        for j = 1:num_shuffle
            rand_ind = datasample(1:size(fr_train_bin_shuffle,1),size(fr_train_bin_shuffle,1));
            fr_temp = fr_train_bin_shuffle(rand_ind,:);
            lr_temp = lr_train_bin_shuffle(rand_ind,:);
            [BS_fr,BS_lr] = get_trains_for_xcross (fr_temp, lr_temp);
            [r_temp,lags] = xcorr(BS_fr,BS_lr,'normalized');
            r_BS_all(j,:) = r_temp([lags>=lag_window(1)/bin_size & lags<=lag_window(2)/bin_size]);
            r_BS(j) = unique(max(r_BS_all(j,:)));
        end

        % xcorr is significant if the 1st percentile of the bootstrapped r is larger than the median shuffied r
        sig_neuron = prctile(r_BS,1) > median(r_shuff);        
        if contains(session_name{n},'right') && i == 1 % flip left/right for right sessions to make everything ipsi vs contra
            selectivity_temp(3) = sig_neuron;
        elseif contains(session_name{n},'right') && i == 3
            selectivity_temp(1) = sig_neuron;
        else
            selectivity_temp(i) = sig_neuron;
        end

        % store the distributions of lags, xcorr coeffs and shuffled coeffs
        lag = lag_window_binned([r == max(r)]);
        if ~isempty(lag)
            lag_aiming(n,i) = lag([abs(lag) == min(abs(lag))]);
            r_aiming(n,i) = r_actual;
            r_shuff_aiming(n,i) = median(r_shuff);
        end
    end
    selectivity_aiming(n,:) = selectivity_temp;
    n
end
if plot_individual == 1
    plot_xcorr_by_neuron(plot_fr,plot_lr,plot_xcorr,r_shuff_aiming);
end
end

function [fr_train,lr_train] = get_trains_for_xcross(fr_maxtrix,lr_matrix)
% average firing rate and lick rate across trials
fr_train_temp = nan(1,size(fr_maxtrix,2));
lr_train_temp = nan(1,size(fr_maxtrix,2));
for k = 1:size(fr_maxtrix,2)
    if sum(~isnan(fr_maxtrix(:,k))) >= 10
        fr_train_temp(1,k) = mean(fr_maxtrix(:,k),'omitnan')*1000;
        lr_train_temp(1,k) = mean(lr_matrix(:,k),'omitnan');
    end
end
% perform mean substraction to prepare for xcorr
fr_train = fr_train_temp-mean(fr_train_temp,'omitnan');
fr_train = fr_train(~isnan(fr_train));
lr_train = lr_train_temp-mean(lr_train_temp,'omitnan');
lr_train = lr_train(~isnan(lr_train));
end

function [] = plot_xcorr_by_neuron(plot_fr,plot_lr,plot_xcorr,r_shuff_aiming)
colors = {'b','r','g'};
for i = 1:size(plot_fr,1)
    if sum(isnan(plot_fr(i,1,:))) ~= size(plot_fr,3)
        figure
        hold on      
        for j = 1:size(plot_fr,2) 
            subplot(3,1,1)
            plot(smoothdata(squeeze(plot_fr(i,j,:)),'movmean',[3 0]),colors{j})
            title(strcat(num2str(i)))
            xlim([1 35])
            hold on
            subplot(3,1,2)
            plot(smoothdata(squeeze(plot_lr(i,j,:)),'movmean',[3 0])*(0.06^3),colors{j})
            xlim([1 35])
            ylim([-25 25])
            hold on
            subplot(3,1,3)
            temp_xcorr = squeeze(plot_xcorr(i,j,:));
            plot(temp_xcorr,colors{j})
            hold on
            if ~isnan(sum(r_shuff_aiming(i,j)))
                xline(find(temp_xcorr==max(temp_xcorr)),colors{j})
                yline(r_shuff_aiming(i,j),colors{j})
                hold on
            end
            xlim([1 21])
            xticks([1 11 21])
            xticklabels({'-100','0','100'})
            ylim([-1 1])
        end
    end
end
end

function [] = plot_xcorr_summary(selectivity_aiming,r_aiming,lag_aiming)
% only include neurons with significant xcorr for at least one trial type
selectivity_aiming = selectivity_aiming(:,1:3);
selectivity_all = sum(selectivity_aiming,2,'omitnan');
sig_ind = find(selectivity_all>=1)';

plot_data = [];
for i = sig_ind
   temp_ind = find(r_aiming(i, :) == max(r_aiming(i, :)));
   plot_data = [plot_data,lag_aiming(i,temp_ind)];
end
figure
plot_x = linspace(-105,105,22);
plot_y = histc(plot_data*10,plot_x)/numel(plot_data);
stairs(plot_x,plot_y,'Color','k')
xlim([-105 105])
ylim([0 0.2])
xlabel('Lag(ms)')
ylabel('Probability')

plot_data_r = [];
plot_data_r_shuff = [];
for i = sig_ind
   plot_data_r = [plot_data_r,max(r_aiming(i, :))];
end
figure
plot_x = linspace(-1,1,21);
plot_y = histc(plot_data_r,plot_x)/numel(plot_data_r);
stairs(plot_x,plot_y,'Color','k')
xlim([0 1.05])
ylim([0 0.5])
xlabel('Max x-corr coef')
ylabel('Probability')
end
%
function [venn_data] = get_venn_data(selectivity_aiming)
for i = 1:size(selectivity_aiming,1)
    if selectivity_aiming(i,1) == 1 && selectivity_aiming(i,2) ~= 1 && selectivity_aiming(i,3) ~= 1
        selectivity_aiming(i,4) = 1;
    elseif selectivity_aiming(i,1) ~= 1 && selectivity_aiming(i,2) == 1 && selectivity_aiming(i,3) ~= 1
        selectivity_aiming(i,4) = 2;
    elseif selectivity_aiming(i,1) ~= 1 && selectivity_aiming(i,2) ~= 1 && selectivity_aiming(i,3) == 1
        selectivity_aiming(i,4) = 3;
    elseif selectivity_aiming(i,1) == 1 && selectivity_aiming(i,2) == 1 && selectivity_aiming(i,3) ~= 1
        selectivity_aiming(i,4) = 4;
    elseif selectivity_aiming(i,1) == 1 && selectivity_aiming(i,2) ~= 1 && selectivity_aiming(i,3) == 1
        selectivity_aiming(i,4) = 5;      
    elseif selectivity_aiming(i,1) ~= 1 && selectivity_aiming(i,2) == 1 && selectivity_aiming(i,3) == 1
        selectivity_aiming(i,4) = 6;
    elseif selectivity_aiming(i,1) == 1 && selectivity_aiming(i,2) == 1 && selectivity_aiming(i,3) == 1
        selectivity_aiming(i,4) = 7;
    else
        selectivity_aiming(i,4) = 0;
    end    
end

for i = 1:7
    venn_data(i) = numel(find(selectivity_aiming(:,4)==i));
end
venn_data = [venn_data, size(selectivity_aiming,1) - sum(venn_data)];
end