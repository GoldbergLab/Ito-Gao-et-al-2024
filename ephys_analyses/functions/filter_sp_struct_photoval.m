function sp_struct_filt = filter_sp_struct_photoval(sp_struct, params)

min_trial_num = params.min_trial_num;
min_fr_thresh = params.min_fr_thresh;

v = 1;
for x = 1:numel(sp_struct)

    % ensure we have enough trials in each condition
    num_laser_trials = sum([sp_struct(x).laser_trial{:}] == 1);
    num_no_laser_trials = sum([sp_struct(x).laser_trial{:}] == 0);
    if num_laser_trials >= min_trial_num && num_no_laser_trials >= min_trial_num

        % calculate mean firing rates across conditions
        sp_times_sec = cellfun(@(x) x/20000, sp_struct(x).sp_times, 'UniformOutput', false);
        mean_firing_rate = mean(cellfun(@(x) sum(x >= -0.3 & x <= 0.3)/0.6, sp_times_sec));

        % ensure mean_firing_rate is > min_fr_thresh
        if mean_firing_rate >= min_fr_thresh

            % filter sp_struct

            % convert spike times from samples to ms
            sp_times = sp_struct(x).sp_times;
            sp_struct(x).sp_times_filt = cellfun(@(z) z/20000*1000, sp_times, 'UniformOutput', false);

            sp_struct_filt(v) = sp_struct(x);
            v = v + 1;
        end
    end
end
end