function [spatial_info_all, spatial_info_shuff, spatial_info_bin, spatial_info_null_all] = skaggs_spatial_info(firing_rate_sig, firing_rate_bin, firing_rate_bin_shuff, num_trials_bin, mean_fr_shuff, mean_fr, lick_num)

spatial_info_all = nan(numel(firing_rate_sig), 1);
spatial_info_null_all = nan(numel(firing_rate_sig), 1);
spatial_info_bin = nan(numel(firing_rate_sig), 1);
spatial_info_shuff = nan(numel(firing_rate_sig), size(firing_rate_bin_shuff, 4), 1);
for cell_num = 1:numel(firing_rate_sig)

    % get total number of trials
    total_trial_num = nansum(num_trials_bin(cell_num, lick_num, :));

    spatial_info = 0;
    spatial_info_all_temp = nan(size(firing_rate_bin, 3), 1);
    for i = 1:size(firing_rate_bin, 3)
        if ~isnan(firing_rate_bin(cell_num, lick_num, i))
            fr_ratio = firing_rate_bin(cell_num, lick_num, i)/mean_fr(cell_num, lick_num);

            % do this to ensure log is valid
            if fr_ratio == 0
                fr_ratio = eps;
            elseif isnan(fr_ratio)
                fr_ratio = eps;
            end

            % calculate occupancy ratio
            occupancy_ratio = num_trials_bin(cell_num, lick_num, i)/total_trial_num;

            % calculate spatial information via Skaggs
            spatial_info_temp = occupancy_ratio*fr_ratio*log2(fr_ratio);
            spatial_info_all_temp(i) = spatial_info_temp;
            spatial_info = spatial_info + spatial_info_temp;
        end
    end
    spatial_info_all(cell_num, 1) = spatial_info;

    % calculate a null for spatial information, where the value of the
    % spatial info at each bin is simply the mean across the respective
    % tertile of the binned distribution. 
    spatial_info_null = 0;
    spatial_info_all_temp = nan(size(firing_rate_bin, 3), 1);
    for i = 1:size(firing_rate_bin, 3)
        if ~isnan(firing_rate_bin(cell_num, lick_num, i))

            bin_ind = size(firing_rate_bin, 3)/3;
            if i <= bin_ind
                firing_rate_bin_temp2 = mean(firing_rate_bin(cell_num, lick_num, 1:5), 'omitnan');
            elseif i > bin_ind && i <= bin_ind*2
                firing_rate_bin_temp2 = mean(firing_rate_bin(cell_num, lick_num, 6:10), 'omitnan');
            elseif i > bin_ind*2
                firing_rate_bin_temp2 = mean(firing_rate_bin(cell_num, lick_num, 11:end), 'omitnan');
            end

            fr_ratio_null = firing_rate_bin_temp2/mean(firing_rate_bin(cell_num, lick_num, :), 'omitnan');
            if fr_ratio_null == 0
                fr_ratio_null = eps;
            elseif isnan(fr_ratio_null)
                fr_ratio_null = eps;
            end

            occupancy_ratio = num_trials_bin(cell_num, lick_num, i)/total_trial_num;

            spatial_info_null_temp = occupancy_ratio*fr_ratio_null*log2(fr_ratio_null);
            spatial_info_all_temp(i) = spatial_info_null_temp;
            spatial_info_null = spatial_info_null + spatial_info_null_temp;
        end
    end
    spatial_info_null_all(cell_num, 1) = spatial_info_null;

    % find the bin with maximal spatial information. If two bins have the
    % exact same spatial info value (very few if any), then take the mean.
    temp = find(spatial_info_all_temp == max(spatial_info_all_temp));
    if numel(temp) == 1
        spatial_info_bin(cell_num, 1) = temp;
    else
        spatial_info_bin(cell_num, 1) = mean(temp);
    end
    
    % calculate the shuffled spatial info for statistical testing
    for i = 1:size(firing_rate_bin_shuff, 4)
        spatial_info_temp_all = 0;
        for j = 1:size(firing_rate_bin_shuff, 3)
            if ~isnan(firing_rate_bin_shuff(cell_num, lick_num, j, i))
                fr_ratio = firing_rate_bin_shuff(cell_num, lick_num, j, i)/mean_fr_shuff(cell_num, lick_num, i);

                if fr_ratio == 0
                    fr_ratio = eps;
                elseif isnan(fr_ratio)
                    fr_ratio = eps;
                end

                occupancy_ratio = num_trials_bin(cell_num, lick_num, j)/total_trial_num;
                spatial_info_temp = occupancy_ratio*fr_ratio*log2(fr_ratio);
                spatial_info_temp_all = spatial_info_temp_all + spatial_info_temp;
            end
        end
        spatial_info_shuff(cell_num, i, 1) = spatial_info_temp_all;
    end
end
end