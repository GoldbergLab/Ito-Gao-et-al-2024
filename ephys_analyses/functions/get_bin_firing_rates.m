function sp_times_bin_align = get_bin_firing_rates(sp_times, align_var, lick_num, align_lick_num, spike_window, bin_size, convert_to_sec)

% create bin windows
window_bin = spike_window(1):bin_size:spike_window(2);

% align sp_times, convert time to sec, get firing rates in window.
sp_times_bin_align = nan(numel(sp_times), numel(window_bin) - 1);
for j = 1:numel(sp_times)

    % if align_var is empty, all spikes are already aligned to cue.  if it
    % is not empty, use align_var to align spikes to.
    if isempty(align_var)
        sp_times_align_temp = sp_times{j}/1000;
    else
        sp_times_align_temp = (sp_times{j} - align_var(lick_num == align_lick_num, j))/1000;
    end

    for k = 1:(numel(window_bin) - 1)  
        % bin number of spikes in window
        bin_temp = sum(sp_times_align_temp >= window_bin(k) & sp_times_align_temp < window_bin(k + 1));

        if convert_to_sec == 1
            % convert to spikes per second
            sp_times_bin_align(j, k) = bin_temp*(1/bin_size);
        else
            sp_times_bin_align(j, k) = bin_temp;
        end
    end
end