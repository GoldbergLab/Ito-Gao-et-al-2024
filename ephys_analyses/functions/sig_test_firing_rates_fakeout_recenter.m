function sig_struct = sig_test_firing_rates_fakeout_recenter(sp_struct, sp_vars_struct, params, lick_num, no_recenter)

if ~exist('no_recenter', 'var')
    no_recenter = false;
end

spike_window = params.spike_window;
fakeout_alpha = params.fakeout_alpha;
fakeout_trial_compare = params.fakeout_trial_compare;
sig_lick_num = params.sig_lick_num;

sp_times_filt = {sp_vars_struct.sp_times_filt};
sp_contact_filt = {sp_vars_struct.sp_contact_filt};
fakeout_trial_filt = {sp_vars_struct.fakeout_trial_filt};

fakeout_p_val = nan(numel(sp_struct), numel(sig_lick_num));
fakeout_sig = nan(numel(sp_struct), numel(sig_lick_num));
fakeout_pairwise_sig = nan(numel(sp_struct), numel(sig_lick_num), 3, 6);
fakeout_mean_firing_rate = nan(numel(sp_struct), numel(sig_lick_num), 3);
fakeout_firing_rate = cell(numel(sp_struct), numel(sig_lick_num), 3);

if no_recenter == false
    recenter_trial_filt = {sp_vars_struct.recenter_trial_filt};
    recenter_p_val = nan(numel(sp_struct), numel(fakeout_trial_compare));
    recenter_sig = nan(numel(sp_struct), numel(fakeout_trial_compare));
    recenter_mean_firing_rate = nan(numel(sp_struct), numel(fakeout_trial_compare), 2);
    recenter_firing_rate = cell(numel(sp_struct), numel(fakeout_trial_compare), 2);
    selectivity = cell(numel(sp_struct), 5);
    selectivity_cent = cell(numel(sp_struct), 5);
else
    selectivity = cell(numel(sp_struct), 3);
    selectivity_cent = cell(numel(sp_struct), 3);
end
firing_rate = cell(numel(sp_struct), 1);
sp_times_align = cell(numel(sp_struct), 1);
for x = 1:numel(sp_struct)      

    % align spike times and get firing rates in window
    [firing_rate_temp, sp_times_align_temp] = get_firing_rates_in_win(sp_times_filt{x}, sp_contact_filt{x}, spike_window, lick_num);
    firing_rate{x} = firing_rate_temp;
    sp_times_align{x} = sp_times_align_temp;

    % Test for significance in firing rates between fakeout and recenter
    % Loop through each lick defined by sig_lick_num (must be included in
    % lick_num), test for significance between fakeout trial types using HSD. 
    for i = 1:numel(sig_lick_num)
        % since we are doing this on recentering, the significance between
        % conditions must be taken only on non-recentering trials for
        % fakeout on licks 4 and beyond
        if sig_lick_num(i) <= 3
            [fp, ~, stats] = anova1(firing_rate_temp(lick_num == sig_lick_num(i), :), fakeout_trial_filt{x}, 'off');
            fakeout_p_val(x, i) = fp;
            if fp < fakeout_alpha
                fakeout_sig(x, i) = 1;
                fakeout_pairwise_sig(x, i, :, :) = multcompare(stats, "CType", "hsd", 'Display', 'off');
            else
                fakeout_sig(x, i) = 0;
            end
        else
            [fp, ~, stats] = anova1(firing_rate_temp(lick_num == sig_lick_num(i), recenter_trial_filt{x} == 0),  fakeout_trial_filt{x}(recenter_trial_filt{x} == 0), 'off');
            fakeout_p_val(x, i) = fp;
            if fp < fakeout_alpha
                fakeout_sig(x, i) = 1;
                fakeout_pairwise_sig(x, i, :, :) = multcompare(stats, "CType", "hsd", 'Display', 'off');
            else
                fakeout_sig(x, i) = 0;
            end 
        end
    end

    % Calculate mean firing rates in each condition
    for i = 1:numel(sig_lick_num)
        for j = 1:size(fakeout_mean_firing_rate, 3)
            if sig_lick_num(i) <= 3
                % x is the neuron, i is the licks in sig_lick_num, j is the
                % fakeout condition (1 is left, 2 center, 3 right).
                fakeout_mean_firing_rate(x, i, j) = mean(firing_rate_temp(lick_num == sig_lick_num(i), fakeout_trial_filt{x} == j));
                fakeout_firing_rate{x, i, j} = firing_rate_temp(lick_num == sig_lick_num(i), fakeout_trial_filt{x} == j);
            else
                % x is the neuron, i is the licks in sig_lick_num, j is the
                % fakeout condition (1 is left, 2 center, 3 right).
                fakeout_mean_firing_rate(x, i, j) = mean(firing_rate_temp(lick_num == sig_lick_num(i), fakeout_trial_filt{x} == j & recenter_trial_filt{x} == 0));
                fakeout_firing_rate{x, i, j} = firing_rate_temp(lick_num == sig_lick_num(i), fakeout_trial_filt{x} == j & recenter_trial_filt{x} == 0);
            end
        end
    end

    if no_recenter == false
        % Test for signficance between recenter conditions for each fakeout trial
        % type using ttest2. For recenter_sig, 0 is not significant, 1 is
        % signficant. If 1 is in first column, it indicates this neuron is 
        % significant for RIGHT nicks (mouse licking left), if 1 is in 2nd column,
        % this neuron is significant for LEFT nicks (mouse licking right).
        recenter_alpha = params.recenter_alpha;
        if sum(lick_num == 4) == 1
            firing_rate_L4 = firing_rate_temp(lick_num == 4, :); 
            for j = 1:numel(fakeout_trial_compare)
                [~, rp] = ttest2(firing_rate_L4(fakeout_trial_filt{x} == fakeout_trial_compare(j) & recenter_trial_filt{x} == 0),...
                    firing_rate_L4(fakeout_trial_filt{x} == fakeout_trial_compare(j) & recenter_trial_filt{x} == 1));
                recenter_p_val(j) = rp;
                if rp < recenter_alpha
                    recenter_sig(x, j) = 1;
                else
                    recenter_sig(x, j) = 0;
                end
    
                % x is the neuron, j is licks in fakeout_trial_compare, 1
                % is no recenter, 2 is recenter.
                recenter_mean_firing_rate(x, j, 1) = mean(firing_rate_L4(fakeout_trial_filt{x} == fakeout_trial_compare(j) & recenter_trial_filt{x} == 0));
                recenter_mean_firing_rate(x, j, 2) = mean(firing_rate_L4(fakeout_trial_filt{x} == fakeout_trial_compare(j) & recenter_trial_filt{x} == 1));
                recenter_firing_rate{x, j, 1} = firing_rate_L4(fakeout_trial_filt{x} == fakeout_trial_compare(j) & recenter_trial_filt{x} == 0);
                recenter_firing_rate{x, j, 2} = firing_rate_L4(fakeout_trial_filt{x} == fakeout_trial_compare(j) & recenter_trial_filt{x} == 1);
            end
        end
    end

    % I could have done this in a number of different ways, there are many
    % 'flavors' of SC neuron responses.  However, all we care about for this
    % analysis is to figure out whether neurons modulated by left or right
    % fakeouts are also modulated by the same nick location during recentering.
    posthoc_alpha = 0.05/3;
    if ~isnan(fakeout_pairwise_sig(x, sig_lick_num == 2, 2, 6))
        % if this neuron is activated by left fakeouts more than right,
        % make it a left selective neuron. Note this does NOT take into
        % account the CENTER response.
        if fakeout_pairwise_sig(x, sig_lick_num == 2, 2, 6) < posthoc_alpha && (fakeout_mean_firing_rate(x, sig_lick_num == 2, 1) > fakeout_mean_firing_rate(x, sig_lick_num == 2, 3))
            selectivity{x, 1} = 1;
        % if this neuron is activated by right fakeouts more than left,
        % then make it a right selective neuron. Note this does NOT take
        % into account the CENTER response.
        elseif fakeout_pairwise_sig(x, sig_lick_num == 2, 2, 6) < posthoc_alpha && (fakeout_mean_firing_rate(x, sig_lick_num == 2, 3) > fakeout_mean_firing_rate(x, sig_lick_num == 2, 1)) 
            selectivity{x, 1} = 3;
        % if this neuron is not significantly different from left vs. right
        elseif fakeout_pairwise_sig(x, sig_lick_num == 2, 2, 6) > posthoc_alpha
            % but this neuron is significant for center vs. left/right,
            % and the firing rate is highest for center, then make it a
            % center selective neuron.
            if fakeout_pairwise_sig(x, sig_lick_num == 2, 1, 6) < posthoc_alpha && fakeout_pairwise_sig(x, sig_lick_num == 2, 3, 6) < posthoc_alpha
                selectivity{x, 1} = 2;
            % below will maybe be only 1 - 2 neurons in the dataset, hard
            % to imagine how this will be the case.  But this will take
            % care of whether left vs. right is not significant, and only
            % left vs. center or left vs. right, but not both, are
            % significant.
            elseif fakeout_pairwise_sig(x, sig_lick_num == 2, 1, 6) < posthoc_alpha || fakeout_pairwise_sig(x, sig_lick_num == 2, 3, 6) < posthoc_alpha
                if (fakeout_mean_firing_rate(x, sig_lick_num == 2, 1) > fakeout_mean_firing_rate(x, sig_lick_num == 2, 2)) && (fakeout_mean_firing_rate(x, sig_lick_num == 2, 1) > fakeout_mean_firing_rate(x, sig_lick_num == 2, 3))
                    selectivity{x, 1} = 1;
                elseif (fakeout_mean_firing_rate(x, sig_lick_num == 2, 3) > fakeout_mean_firing_rate(x, sig_lick_num == 2, 2)) && (fakeout_mean_firing_rate(x, sig_lick_num == 2, 3) > fakeout_mean_firing_rate(x, sig_lick_num == 2, 1))
                    selectivity{x, 1} = 3;
                elseif (fakeout_mean_firing_rate(x, sig_lick_num == 2, 2) > fakeout_mean_firing_rate(x, sig_lick_num == 2, 1)) && (fakeout_mean_firing_rate(x, sig_lick_num == 2, 2) > fakeout_mean_firing_rate(x, sig_lick_num == 2, 3))
                    selectivity{x, 1} = 2;
                end
            else
                selectivity{x, 1} = 0;
            end
        end
    elseif isnan(fakeout_pairwise_sig(x, sig_lick_num == 2, 2, 6))
        selectivity{x, 1} = 0;
    end

    % I could have done this in a number of different ways, there are many
    % 'flavors' of SC neuron responses.  However, all we care about for this
    % analysis is to figure out whether neurons modulated by left or right
    % fakeouts are also modulated by the same nick location during recentering.
    posthoc_alpha = 0.05/3;
    if ~isnan(fakeout_pairwise_sig(x, sig_lick_num == 2, 2, 6))
        % if this neuron is activated by left fakeouts more than right, AND
        % CENTER, then make it a left selective neuron. 
        if fakeout_pairwise_sig(x, sig_lick_num == 2, 2, 6) < posthoc_alpha
            if (fakeout_mean_firing_rate(x, sig_lick_num == 2, 1) > fakeout_mean_firing_rate(x, sig_lick_num == 2, 3))
                if fakeout_pairwise_sig(x, sig_lick_num == 2, 1, 6) < posthoc_alpha && (fakeout_mean_firing_rate(x, sig_lick_num == 2, 1) > fakeout_mean_firing_rate(x, sig_lick_num == 2, 2))
                    selectivity_cent{x, 1} = 1;
                elseif fakeout_pairwise_sig(x, sig_lick_num == 2, 1, 6) < posthoc_alpha && (fakeout_mean_firing_rate(x, sig_lick_num == 2, 2) > fakeout_mean_firing_rate(x, sig_lick_num == 2, 1))
                    selectivity_cent{x, 1} = 2;
                elseif fakeout_pairwise_sig(x, sig_lick_num == 2, 1, 6) > posthoc_alpha
                    selectivity_cent{x, 1} = 4;
                end
            elseif (fakeout_mean_firing_rate(x, sig_lick_num == 2, 3) > fakeout_mean_firing_rate(x, sig_lick_num == 2, 1))
                if fakeout_pairwise_sig(x, sig_lick_num == 2, 3, 6) < posthoc_alpha && (fakeout_mean_firing_rate(x, sig_lick_num == 2, 3) > fakeout_mean_firing_rate(x, sig_lick_num == 2, 2))
                    selectivity_cent{x, 1} = 3;
                elseif fakeout_pairwise_sig(x, sig_lick_num == 2, 3, 6) < posthoc_alpha && (fakeout_mean_firing_rate(x, sig_lick_num == 2, 2) > fakeout_mean_firing_rate(x, sig_lick_num == 2, 3))
                    selectivity_cent{x, 1} = 2;
                elseif fakeout_pairwise_sig(x, sig_lick_num == 2, 3, 6) > posthoc_alpha
                    selectivity_cent{x, 1} = 5;
                end
            end
        % if this neuron is not significantly different from left vs. right
        elseif fakeout_pairwise_sig(x, sig_lick_num == 2, 2, 6) > posthoc_alpha
            % but this neuron is significant for center vs. left/right,
            % and the firing rate is highest for center, then make it a
            % center selective neuron.
            if fakeout_pairwise_sig(x, sig_lick_num == 2, 1, 6) < posthoc_alpha && fakeout_pairwise_sig(x, sig_lick_num == 2, 3, 6) < posthoc_alpha
                selectivity_cent{x, 1} = 2;
            % below will maybe be only 1 - 2 neurons in the dataset, hard
            % to imagine how this will be the case.  But this will take
            % care of whether left vs. right is not significant, and only
            % left vs. center or left vs. right, but not both, are
            % significant.
            elseif fakeout_pairwise_sig(x, sig_lick_num == 2, 1, 6) < posthoc_alpha || fakeout_pairwise_sig(x, sig_lick_num == 2, 3, 6) < posthoc_alpha
                if (fakeout_mean_firing_rate(x, sig_lick_num == 2, 1) > fakeout_mean_firing_rate(x, sig_lick_num == 2, 2)) && (fakeout_mean_firing_rate(x, sig_lick_num == 2, 1) > fakeout_mean_firing_rate(x, sig_lick_num == 2, 3))
                    selectivity_cent{x, 1} = 1;
                elseif (fakeout_mean_firing_rate(x, sig_lick_num == 2, 3) > fakeout_mean_firing_rate(x, sig_lick_num == 2, 2)) && (fakeout_mean_firing_rate(x, sig_lick_num == 2, 3) > fakeout_mean_firing_rate(x, sig_lick_num == 2, 1))
                    selectivity_cent{x, 1} = 3;
                elseif (fakeout_mean_firing_rate(x, sig_lick_num == 2, 2) > fakeout_mean_firing_rate(x, sig_lick_num == 2, 1)) && (fakeout_mean_firing_rate(x, sig_lick_num == 2, 2) > fakeout_mean_firing_rate(x, sig_lick_num == 2, 3))
                    selectivity_cent{x, 1} = 2;
                end
            else
                selectivity_cent{x, 1} = 0;
            end
        end
    elseif isnan(fakeout_pairwise_sig(x, sig_lick_num == 2, 2, 6))
        selectivity_cent{x, 1} = 0;
    end

    if no_recenter == false     
        % and this neuron fires more for LEFT nicks on recenter (mouse licking
        % right), then assign it a 1. 
        if recenter_sig(x, 2) == 1 && (recenter_mean_firing_rate(x, 2, 2) > recenter_mean_firing_rate(x, 2, 1))
            selectivity{x, 2} = 1;
            selectivity_cent{x, 2} = 1;
        elseif recenter_sig(x, 2) == 1 && (recenter_mean_firing_rate(x, 2, 2) < recenter_mean_firing_rate(x, 2, 1))
            selectivity{x, 2} = -1;
            selectivity_cent{x, 2} = -1;
        % if this neuron is not modulated by LEFT nicks on recentering, than assign it 0.
        elseif recenter_sig(x, 2) == 0
            selectivity{x, 2} = 0;
            selectivity_cent{x, 2} = 0;
        end 
    
        % and this neuron fires more for RIGHT nicks on recenter (mouse licking
        % left), then assign it a 2.
        if recenter_sig(x, 1) == 1 && (recenter_mean_firing_rate(x, 1, 2) > recenter_mean_firing_rate(x, 1, 1))
            selectivity{x, 3} = 3;  
            selectivity_cent{x, 3} = 3;
        elseif recenter_sig(x, 1) == 1 && (recenter_mean_firing_rate(x, 1, 2) < recenter_mean_firing_rate(x, 1, 1))
            selectivity{x, 3} = -3;
            selectivity_cent{x, 3} = -3;
        % if this neuron is not modulated by recentering, than assign it 0.
        elseif recenter_sig(x, 1) == 0
            selectivity{x, 3} = 0;
            selectivity_cent{x, 3} = 0;
        end      

        selectivity{x, 4} = sp_struct(x).sp_cluster_id;
        selectivity{x, 5} = sp_struct(x).session_id;
        selectivity_cent{x, 4} = sp_struct(x).sp_cluster_id;
        selectivity_cent{x, 5} = sp_struct(x).session_id;
    else
        selectivity_cent{x, 2} = sp_struct(x).sp_cluster_id;
        selectivity_cent{x, 3} = sp_struct(x).session_id;
        selectivity{x, 2} = sp_struct(x).sp_cluster_id;
        selectivity{x, 3} = sp_struct(x).session_id;
    end
end

sig_struct = struct('fakeout_p_val', fakeout_p_val);
[sig_struct.fakeout_sig] = fakeout_sig;
[sig_struct.fakeout_pairwise_sig] = fakeout_pairwise_sig;
[sig_struct.fakeout_mean_firing_rate] = fakeout_mean_firing_rate;
[sig_struct.fakeout_firing_rate] = fakeout_firing_rate;
if no_recenter == false
    [sig_struct.recenter_p_val] = recenter_p_val;
    [sig_struct.recenter_sig] = recenter_sig;
    [sig_struct.recenter_mean_firing_rate] = recenter_mean_firing_rate;
    [sig_struct.recenter_firing_rate] = recenter_firing_rate;
end
[sig_struct.selectivity] = selectivity;
[sig_struct.selectivity_cent] = selectivity_cent;
[sig_struct.firing_rate] = firing_rate;
[sig_struct.sp_times_align] = sp_times_align;

end