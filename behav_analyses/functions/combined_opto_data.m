% this code assumes data_set_all is organized as alternating left and right
% unilateral inactivation sessions
function angle_sum = combined_opto_data(data_set_all,lick_num,recentering,ylabel_input,lick_index_contact_num,L2contact_only,mm_pix)

struct_combined = cell(1,numel(data_set_all));
for i = 1:numel(data_set_all)
    load(data_set_all{i})
    struct_combined{i} = t_stats_combined;
end

% load(data_set_one)
% t_stats_combined_one = t_stats_combined;
% load(data_set_two)
% t_stats_combined_two = t_stats_combined;
% struct_combined = {t_stats_combined_one,t_stats_combined_two};
angle_sum = {};

for p = 1:numel(struct_combined)
    prot_angle = cell(numel(struct_combined{p}), 2, 3, numel(lick_num));
    for q = 1:numel(struct_combined{p})
        t_stats = struct_combined{p}(q).data;

        % laser only labels licks when laser was on - need to create another
        % list that labels all licks on a trial when laser was on
        laser_trial_num = unique([t_stats([t_stats.laser_2D] == 1).trial_num]);    
        laser_trial = ismember([t_stats.trial_num], laser_trial_num);

        prot_angle_all_session = zeros(numel(t_stats), 1);
        for i = 1:numel(t_stats)
            % calculate different kinmetiac parameters, variables in this script are all named as angles but they can be substituted with other kinematic parameters
            if strcmp(ylabel_input,'protAngle')
                hypotenuse = sqrt((t_stats(i).tip_xf(t_stats(i).prot_ind)^2) + (t_stats(i).tip_yf(t_stats(i).prot_ind)^2));
                prot_angle_all_session(i) = asind(t_stats(i).tip_xf(t_stats(i).prot_ind)/hypotenuse);
            elseif strcmp(ylabel_input,'latDisp')
                prot_angle_all_session(i) = t_stats(i).tip_xf(t_stats(i).prot_ind)*mm_pix;
            elseif strcmp(ylabel_input,'duration')
                prot_angle_all_session(i) = t_stats(i).dur;
            elseif strcmp(ylabel_input,'pathlength')
                prot_angle_all_session(i) = t_stats(i).pathlength_3D_t*mm_pix;
            elseif strcmp(ylabel_input,'maxSpeed')
                prot_angle_all_session(i) = max(t_stats(i).magspeed_t)*mm_pix*1000;
            elseif strcmp(ylabel_input,'accPeakNum')
                prot_angle_all_session(i) = numel(t_stats(i).accel_peaks_tot_t);
            elseif strcmp(ylabel_input,'InterLickInterval')
                if i>1
                    prot_angle_all_session(i) = t_stats(i).time_rel_cue - t_stats(i-1).time_rel_cue;
                else
                    prot_angle_all_session(i) = NaN;
                end
            elseif strcmp(ylabel_input,'contactON_lickOFF')
                if i < numel(t_stats) && ~isnan(t_stats(i).spout_contact)
                    prot_angle_all_session(i) = t_stats(i).time_rel_cue + t_stats(i).dur - t_stats(i).spout_contact;
                end
            elseif strcmp(ylabel_input,'contactON_protOFF')
                if i>1 && ~isnan(t_stats(i-1).spout_contact)
                    prot_angle_all_session(i) = t_stats(i).time_rel_cue - t_stats(i-1).spout_contact;
                else
                    prot_angle_all_session(i) = NaN;
                end
            elseif strcmp(ylabel_input,'contactSite')
                if ~isnan(t_stats(i).spout_contact) && numel(t_stats(i).contact_centroid) > 1
                    contact_x_median = median(t_stats(i).contact_centroid(:,1)) + unique(t_stats(i).tip_xf-t_stats(i).tip_x);
                    contact_y_median = median(t_stats(i).contact_centroid(:,2)) + unique(t_stats(i).tip_yf-t_stats(i).tip_y);
                    x2=t_stats(i).tip_xf(t_stats(i).prot_ind);
                    y2=t_stats(i).tip_yf(t_stats(i).prot_ind);
                    x1=contact_x_median;
                    y1=contact_y_median;
                    prot_angle_all_session(i) = atan2d(x1*y2-y1*x2,x1*x2+y1*y2); % angle between two vectors
                else
                    prot_angle_all_session(i) = NaN;
                end
            end
        end

        % cell arrays organized by condition (1 - control, 2 - laser), spout
        % position (1/2/3, left/center/right), and lick number (1, 2, 3, etc.)
        prot_angle_temp = cell(2, max([t_stats.fakeout_trial]), numel(lick_num));

        for j = 1:size(prot_angle_temp, 2)
            for k = 1:size(prot_angle_temp, 3)
                if lick_index_contact_num == 1
                    lick_index_contact = [t_stats.lick_index_contact1];
                    % filter trials with L2 contact, if requested
                    if L2contact_only == 1
                        trial_num_L2contact = [t_stats(~isnan([t_stats.spout_contact]) & [t_stats.lick_index_contact1] == 2).trial_num];
                    elseif L2contact_only == 0
                        trial_num_L2contact = [t_stats.trial_num];
                    end

                    % filter trials
                    trial_num_L1_double_tap = [t_stats(isnan([t_stats.spout_contact2]) & [t_stats.lick_index_contact1] == 1).trial_num];
                    trial_num_L2_double_tap = [t_stats(isnan([t_stats.spout_contact2]) & [t_stats.lick_index_contact1] == 2).trial_num];
                    trial_num_L3_double_tap = [t_stats(isnan([t_stats.spout_contact2]) & [t_stats.lick_index_contact1] == 3).trial_num];
                    trial_num_L4_double_tap = [t_stats(isnan([t_stats.spout_contact2]) & [t_stats.lick_index_contact1] == 4).trial_num];
                    trial_num_L3contact = [t_stats(~isnan([t_stats.spout_contact]) & [t_stats.lick_index_contact1] == 3).trial_num];
                    trial_num_L4contact = [t_stats(~isnan([t_stats.spout_contact]) & [t_stats.lick_index_contact1] == 4).trial_num];
                    trial_num_half = [t_stats.trial_num] <= (max([t_stats.trial_num]));
             
                    % calculate the values for histogram & store values in appropriate variables - angle_temp is organized by condition (1 - control, 2 - laser)
                    angle_temp = cell(2, 1);
                    for m = 1:numel(angle_temp)
                        if lick_num(k) < 3 % for L1 and L2
                            if isfield(t_stats, 'spout_contact2')
                                if lick_num(k) == 1
                                    angle_temp(m) = {prot_angle_all_session([t_stats.fakeout_trial] == j & lick_index_contact == lick_num(k) & laser_trial == (m-1) & trial_num_half)};
                                else
                                    angle_temp(m) = {prot_angle_all_session([t_stats.fakeout_trial] == j & ismember([t_stats.trial_num], trial_num_L2contact) & lick_index_contact == lick_num(k) & laser_trial == (m-1) & ismember([t_stats.trial_num], trial_num_L2_double_tap) & ismember([t_stats.trial_num], trial_num_L1_double_tap) & trial_num_half)};
                                end
                            else
                                disp('Warning - no spout_contact2 field in t_stats!')
                                angle_temp(m) = {prot_angle_all_session([t_stats.fakeout_trial] == j & ismember([t_stats.trial_num], trial_num_L2contact) & lick_index_contact == lick_num(k) & laser_trial == (m-1) & trial_num_half)};
                            end
                        elseif recentering == 1 && lick_num(k) == 3 % for L3 in recentering analyses
                            if j == 2
                                angle_temp(m) = {prot_angle_all_session([t_stats.fakeout_trial] == j & ismember([t_stats.trial_num], trial_num_L2contact) & ismember([t_stats.trial_num], trial_num_L3contact) & lick_index_contact == lick_num(k) & ismember([t_stats.trial_num], trial_num_L2_double_tap) & ismember([t_stats.trial_num], trial_num_L1_double_tap) & ismember([t_stats.trial_num], trial_num_L3_double_tap) & trial_num_half )};
                            else    
                                angle_temp(m) = {prot_angle_all_session([t_stats.fakeout_trial] == j & [t_stats.recenter_trial] == (m-1) & ismember([t_stats.trial_num], trial_num_L2contact) & ismember([t_stats.trial_num], trial_num_L3contact) & lick_index_contact == lick_num(k) & ismember([t_stats.trial_num], trial_num_L2_double_tap) & ismember([t_stats.trial_num], trial_num_L1_double_tap) & ismember([t_stats.trial_num], trial_num_L3_double_tap) & trial_num_half )};
                            end
                        elseif recentering == 1 && lick_num(k) > 3 % for L4 and beyond in recentering analyses
                            if j == 2
                                angle_temp(m) = {prot_angle_all_session([t_stats.fakeout_trial] == j & ismember([t_stats.trial_num], trial_num_L2contact) & ismember([t_stats.trial_num], trial_num_L3contact) & ismember([t_stats.trial_num], trial_num_L4contact) & lick_index_contact == lick_num(k) & ismember([t_stats.trial_num], trial_num_L2_double_tap) & ismember([t_stats.trial_num], trial_num_L1_double_tap) & ismember([t_stats.trial_num], trial_num_L3_double_tap) & ismember([t_stats.trial_num], trial_num_L4_double_tap) & trial_num_half )};
                            else    
                                angle_temp(m) = {prot_angle_all_session([t_stats.fakeout_trial] == j & [t_stats.recenter_trial] == (m-1) & ismember([t_stats.trial_num], trial_num_L2contact) & ismember([t_stats.trial_num], trial_num_L3contact) & ismember([t_stats.trial_num], trial_num_L4contact) & lick_index_contact == lick_num(k) & ismember([t_stats.trial_num], trial_num_L2_double_tap) & ismember([t_stats.trial_num], trial_num_L1_double_tap) & ismember([t_stats.trial_num], trial_num_L3_double_tap) & ismember([t_stats.trial_num], trial_num_L4_double_tap) & trial_num_half )};
                            end
                        elseif recentering == 0 && lick_num(k) == 3 % for L3 in non-recentering analyses
                            angle_temp(m) = {prot_angle_all_session([t_stats.fakeout_trial] == j & ismember([t_stats.trial_num], trial_num_L2contact) & lick_index_contact == lick_num(k) & [t_stats.laser_2D] == (m-1) & laser_trial == (m-1) & ismember([t_stats.trial_num], trial_num_L2_double_tap) & ismember([t_stats.trial_num], trial_num_L1_double_tap) & trial_num_half )};
                        elseif recentering == 0 && lick_num(k) > 3
                            angle_temp(m) = {prot_angle_all_session([t_stats.fakeout_trial] == j & ismember([t_stats.trial_num], trial_num_L2contact) & lick_index_contact == lick_num(k) & laser_trial == (m-1) & ismember([t_stats.trial_num], trial_num_L2_double_tap) & ismember([t_stats.trial_num], trial_num_L1_double_tap) & trial_num_half)};
                        end

                        angle_temp{m} = angle_temp{m}(~isnan(angle_temp{m})); % To filter out NaNs. This pretty much only affects contact site analyses.
                        prot_angle_temp{m, j, k} = angle_temp{m};
                    end
                end
            end
        end
        prot_angle(q, :, :, :) = prot_angle_temp;
    end

    med_prot_angle = cellfun(@median, prot_angle);
    angle_sum_temp = cell(1,size(med_prot_angle,3));
    for i = 1:size(med_prot_angle,3)
        angle_sum_temp{i} = squeeze(prot_angle(:,:,i));
    end  
    angle_sum = vertcat(angle_sum,angle_sum_temp);
end
end