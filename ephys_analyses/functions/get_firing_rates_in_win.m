function [firing_rate, sp_times_align] = get_firing_rates_in_win(sp_times, align_var, window, lick_num)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function aligns spike times to align_var, takes the spike counts in
% a window relative to align_var, for specific lick numbers.
%
% align_var is a variable from t_stats with time relative to cue in ms, 
% such as sp_contact, protrusion, etc.
%
% window is an array whose entries specify two time points in which to take
% spike counts, in seconds. ex., [0 0.1] will take a 100ms snippet after
% the event you are aligning spikes to.
%
% lick_num is an array whose entries specify the licks to align spikes to
% relative to align_var. ex., if align_var is sp_contact, lick num == [2 3
% 4] will take spike counts relative to 
%
% Author:  Brendan Ito
% Email:   bsi8@cornell.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% align sp_times, convert time to sec, get firing rates in window.
for j = 1:numel(sp_times)
    for k = 1:numel(lick_num)
        sp_times_align_temp = (sp_times{j} - align_var(lick_num(k), j))/1000;
        sp_times_align(k, j) = {sp_times_align_temp};

        if ~isempty(window)
            sp_times_in_win = sp_times_align_temp(sp_times_align_temp <= window(2) & sp_times_align_temp > window(1));
            if ~isempty(sp_times_in_win)
                sp_times_align_win(k, j) = {sp_times_in_win};
                firing_rate(k, j) = numel(sp_times_in_win)/(window(2) - window(1));
            elseif isempty(sp_times_in_win)
                sp_times_align_win(k, j) = {0};
                firing_rate(k, j) = 0;
            end
        else
            firing_rate = [];
        end
    end
end