function [sp_vars_struct, sp_struct_filt] = filter_sp_neurons_and_vars_recenter(sp_struct, min_trial_num, min_fr_thresh, lick_num, params, no_recenter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function takes an sp_struct file for sessions and outputs a filtered 
% sp_sturct file (sp_struct_filt), which filters out neurons depending on 
% whether they have more than min_trial_num in each condition, and 
% whether activity of neurons during trial period are > min_fr_thresh.
%
% Author:  Brendan Ito
% Email:   bsi8@cornell.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('no_recenter', 'var')
    no_recenter = false;
end

lick_exist_lick_num = params.lick_exist_lick_num;
sp_contact_lick_num = params.sp_contact_lick_num;
sp_contact2_absent_lick_num = params.sp_contact2_absent_lick_num;
contact_centroid_lick_num = params.contact_centroid_lick_num;

v = 1;
for x = 1:numel(sp_struct)

    lick_index_contact1 = sp_struct(x).lick_index_contact1;

    % filter trials in sp_struct 
    filter_ind = filter_trials_sp_struct(sp_struct, x, lick_index_contact1, lick_exist_lick_num,...
        sp_contact_lick_num, sp_contact2_absent_lick_num, contact_centroid_lick_num);

    % after pre-filtering, check if each condition has min num of trials
    if no_recenter == false
        % after pre-filtering, check if each recentering condition has the minimum number of trials
        fakeout_temp = sp_struct(x).fakeout_trial;
        fakeout_filt_temp = cellfun(@(x) x(1), fakeout_temp(filter_ind));
        recenter_temp = sp_struct(x).recenter_trial;
        recenter_filt_temp = cellfun(@(x) x(1), recenter_temp(filter_ind));
        for j = 1:max(fakeout_filt_temp)
            if j == 1 || j == 3
                for k = 1:2
                    trial_min_filt_temp(j, k) = sum(fakeout_filt_temp == j & recenter_filt_temp == k - 1) > min_trial_num;
                end
            elseif j == 2
                trial_min_filt_temp(j, 1:2) = sum(fakeout_filt_temp == j) >= min_trial_num;
            end
        end
    elseif no_recenter == true
        fakeout_temp = sp_struct(x).fakeout_trial;
        fakeout_filt_temp = cellfun(@(x) x(1), fakeout_temp(filter_ind));
        recenter_temp = sp_struct(x).recenter_trial;
        recenter_filt_temp = cellfun(@(x) x(1), recenter_temp(filter_ind));
        for j = 1:max(fakeout_filt_temp)
            trial_min_filt_temp(j, 1:2) = sum(fakeout_filt_temp == j) >= min_trial_num;
        end
    end

    % calculate mean firing rates across conditions
    sp_times_sec = cellfun(@(x) x/20000, sp_struct(x).sp_times, 'UniformOutput', false);
    mean_firing_rate = mean(cellfun(@(x) sum(x > 0 & x < 1.3)/1.3, sp_times_sec));

    % if all conditions pass...
    if sum(trial_min_filt_temp, 'all') == 6 && mean_firing_rate > min_fr_thresh

        % filter sp_struct
        sp_struct_filt(v) = sp_struct(x);

        % filter relevant sp_struct variables
        fakeout_trial_filt{v} = fakeout_filt_temp;
        recenter_trial_filt{v} = recenter_filt_temp;
        lick_index_contact1_filt{v} = lick_index_contact1(filter_ind);
    
        % convert spike times from samples to ms
        sp_times = sp_struct(x).sp_times;
        sp_times_filt{v} = cellfun(@(z) z/20000*1000, sp_times(filter_ind), 'UniformOutput', false);
                
        % extract kinematic variables to be filtered
        contact_centroid = sp_struct(x).contact_centroid;
        protrusion_onset = sp_struct(x).protrusion;
        protrusion_offset = sp_struct(x).protrusion_offset;
        protrusion_offset_cue = cellfun(@(x, y) x + y, protrusion_onset, protrusion_offset, 'UniformOutput', false);
        sp_contact = sp_struct(x).sp_contact;
        tip_xf_prot = sp_struct(x).tip_xf_prot;
        tip_yf_prot = sp_struct(x).tip_yf_prot;
        tip_zf_prot = sp_struct(x).tip_zf_prot;
        tip_xf_cont = sp_struct(x).tip_xf_cont;
        tip_yf_cont = sp_struct(x).tip_yf_cont;
        tip_zf_cont = sp_struct(x).tip_zf_cont;
        centroid_xf_cont = sp_struct(x).centroid_xf_cont;
        centroid_yf_cont = sp_struct(x).centroid_yf_cont;
        centroid_zf_cont = sp_struct(x).centroid_zf_cont;
        dur = sp_struct(x).dur;
        prot_off = sp_struct(x).protrusion_offset;
        peak_prot_speed = sp_struct(x).peak_prot_speed;
        volume = sp_struct(x).volume;
        if isfield(sp_struct, 'contact_dur')
            contact_dur = sp_struct(x).contact_dur;
        end
        
        % filter relevant kinematic variables
        for i = 1:numel(lick_num)
            prot_off_cue_filt{v}(i, :) = cellfun(@(z, y) y(z == lick_num(i)), lick_index_contact1_filt{v}, protrusion_offset_cue(filter_ind));
            sp_contact_filt{v}(i, :) = cellfun(@(z, y) y(z == lick_num(i)), lick_index_contact1_filt{v}, sp_contact(filter_ind));
    
            contact_centroid_filt{v}(i, :) = cellfun(@(z, y) y(z == lick_num(i)), lick_index_contact1_filt{v}, contact_centroid(filter_ind));
            tip_xf_prot_filt{v}(i, :) = cellfun(@(z, y) y(z == lick_num(i)), lick_index_contact1_filt{v}, tip_xf_prot(filter_ind));
            tip_yf_prot_filt{v}(i, :) = cellfun(@(z, y) y(z == lick_num(i)), lick_index_contact1_filt{v}, tip_yf_prot(filter_ind));
            tip_zf_prot_filt{v}(i, :) = cellfun(@(z, y) y(z == lick_num(i)), lick_index_contact1_filt{v}, tip_zf_prot(filter_ind));
            tip_xf_cont_filt{v}(i, :) = cellfun(@(z, y) y(z == lick_num(i)), lick_index_contact1_filt{v}, tip_xf_cont(filter_ind));
            tip_yf_cont_filt{v}(i, :) = cellfun(@(z, y) y(z == lick_num(i)), lick_index_contact1_filt{v}, tip_yf_cont(filter_ind));
            tip_zf_cont_filt{v}(i, :) = cellfun(@(z, y) y(z == lick_num(i)), lick_index_contact1_filt{v}, tip_zf_cont(filter_ind));
            centroid_xf_cont_filt{v}(i, :)  = cellfun(@(z, y) y(z == lick_num(i)), lick_index_contact1_filt{v}, centroid_xf_cont(filter_ind));
            centroid_yf_cont_filt{v}(i, :)  = cellfun(@(z, y) y(z == lick_num(i)), lick_index_contact1_filt{v}, centroid_yf_cont(filter_ind));
            centroid_zf_cont_filt{v}(i, :)  = cellfun(@(z, y) y(z == lick_num(i)), lick_index_contact1_filt{v}, centroid_zf_cont(filter_ind));
            dur_filt{v}(i, :) = cellfun(@(z, y) y(z == lick_num(i)), lick_index_contact1_filt{v}, dur(filter_ind));
            prot_off_filt{v}(i, :) = cellfun(@(z, y) y(z == lick_num(i)), lick_index_contact1_filt{v}, prot_off(filter_ind));
            peak_prot_speed_filt{v}(i, :) = cellfun(@(z, y) y(z == lick_num(i)), lick_index_contact1_filt{v}, peak_prot_speed(filter_ind));
            volume_filt{v}(i, :) = cellfun(@(z, y) y(z == lick_num(i)), lick_index_contact1_filt{v}, volume(filter_ind));
            prot_onset_filt{v}(i, :) = cellfun(@(z, y) y(z == lick_num(i)), lick_index_contact1_filt{v}, protrusion_onset(filter_ind));
            if isfield(sp_struct, 'contact_dur')
                contact_dur_filt{v}(i, :) = cellfun(@(z, y) y(z == lick_num(i)), lick_index_contact1_filt{v}, contact_dur(filter_ind));
            end
        end

        v = v + 1;

    end
end

% assign vars to fields in sp_vars_struct
sp_vars_struct = struct('fakeout_trial_filt', fakeout_trial_filt);
[sp_vars_struct.recenter_trial_filt] = deal(recenter_trial_filt{:});
[sp_vars_struct.lick_index_contact1_filt] = deal(lick_index_contact1_filt{:});
[sp_vars_struct.sp_times_filt] = deal(sp_times_filt{:});
[sp_vars_struct.prot_off_cue_filt] = deal(prot_off_cue_filt{:});
[sp_vars_struct.sp_contact_filt] = deal(sp_contact_filt{:});
[sp_vars_struct.contact_centroid_filt] = deal(contact_centroid_filt{:});
[sp_vars_struct.tip_xf_prot_filt] = deal(tip_xf_prot_filt{:});
[sp_vars_struct.tip_yf_prot_filt] = deal(tip_yf_prot_filt{:});
[sp_vars_struct.tip_zf_prot_filt] = deal(tip_zf_prot_filt{:});
[sp_vars_struct.tip_xf_cont_filt] = deal(tip_xf_cont_filt{:});
[sp_vars_struct.tip_yf_cont_filt] = deal(tip_yf_cont_filt{:});
[sp_vars_struct.tip_zf_cont_filt] = deal(tip_zf_cont_filt{:});
[sp_vars_struct.centroid_xf_cont_filt] = deal(centroid_xf_cont_filt{:});
[sp_vars_struct.centroid_yf_cont_filt] = deal(centroid_yf_cont_filt{:});
[sp_vars_struct.centroid_zf_cont_filt] = deal(centroid_zf_cont_filt{:});
[sp_vars_struct.dur_filt] = deal(dur_filt{:});
[sp_vars_struct.prot_off_filt] = deal(prot_off_filt{:});
[sp_vars_struct.peak_prot_speed_filt] = deal(peak_prot_speed_filt{:});
[sp_vars_struct.volume_filt] = deal(volume_filt{:});
[sp_vars_struct.prot_onset_filt] = deal(prot_onset_filt{:});
if isfield(sp_struct, 'contact_dur')
    [sp_vars_struct.contact_dur_filt] = deal(contact_dur_filt{:});
end

end