function [plot_data_summary,plot_data_indiv,shuffle_test_lesion_tailOne,shuffle_test_lesion_tailTwo,shuffle_test_position_tailOne,shuffle_test_position_tailTwo,shuffle_test_lesion_combined_tailOne,shuffle_test_lesion_combined_tailTwo,plot_data_summary_combined,plot_data_indiv_combined] = plot_xsession_comparison(data_set_one,data_set_two,lick_num,true_L1,ylabel_input,ylim_val,min_max_angle,plot_individual,lick_index_contact_num,L2contact_only,mm_pix,bin_size, fig_name)

load(data_set_one)
animal_t_stats_paths_preLesion = t_stats_combined;
load(data_set_two)
animal_t_stats_paths_postLesion = t_stats_combined;

prot_angle = cell(numel(animal_t_stats_paths_preLesion), 2, 3, numel(lick_num));
med_prot_angle_all = cell(numel(animal_t_stats_paths_preLesion), 2, 3, numel(lick_num));
prot_angle_all_lesion = cell(numel(animal_t_stats_paths_preLesion), 3, numel(lick_num));
prot_angle_all_position = cell(numel(animal_t_stats_paths_preLesion), 2, numel(lick_num));

for p = 1:numel(animal_t_stats_paths_preLesion)
    t_stats = animal_t_stats_paths_preLesion(p).data;
    for i = 1:numel(t_stats)
        t_stats(i).lesion_trial = 0;
    end
    t_stats_preLesion = t_stats;
    t_stats = animal_t_stats_paths_postLesion(p).data;
    for i = 1:numel(t_stats)
        t_stats(i).trial_num = t_stats(i).trial_num + max([t_stats_preLesion.trial_num]);
        t_stats(i).lesion_trial = 1;
    end    
    t_stats = [t_stats_preLesion,t_stats];
       
    prot_angle_all_session = zeros(numel(t_stats), 1);
    for i = 1:numel(t_stats)
        % calculate different kinmetiac parameters, variables in this script are all named as angles but they can be substituted with other kinematic parameters
        
        if isnan(t_stats(i).dur)
           prot_angle_all_session(i) = NaN;
           t_stats(i).lick_index = NaN;
        else
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
            elseif strcmp(ylabel_input,'reaction_time')
                prot_angle_all_session(i) = t_stats(i).time_rel_cue;
            end
        end
    end

    % cell arrays organized by condition (1 - pre-lesion, 2 - post-lesion), spout
    % position (1/2/3, left/center/right), and lick number (1, 2, 3, etc.)
    prot_angle_bin = cell(2, max([t_stats.fakeout_trial]), numel(lick_num));
    prot_angle_temp = cell(2, max([t_stats.fakeout_trial]), numel(lick_num));
    x_bin = linspace(min_max_angle(1), min_max_angle(2), (sum(abs(min_max_angle))+bin_size)/bin_size);
    
    for j = 1:size(prot_angle_bin, 2)
        for k = 1:size(prot_angle_bin, 3)
            if lick_index_contact_num == 1
                lick_index_contact = [t_stats.lick_index_contact1];
                % filter trials with L2 contact, if requested
                if L2contact_only == 0
                    trial_num_L2contact = [t_stats(~isnan([t_stats.spout_contact]) & [t_stats.lick_index_contact1] == 2).trial_num];
                else
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
                
                % angle_bin_temp organized by condition (1 - pre-lesion, 2 - post-lesion), number of bins)
                angle_bin_temp = zeros(2, numel(x_bin));                
                % calculate the values for histogram & store values in appropriate variables - angle_temp is organized by condition (1 - pre-lesion, 2 - post-lesion)
                angle_temp = cell(2, 1);
                for m = 1:numel(angle_temp)
                    if lick_num(k) < 3 % for L1 and L2
                        if isfield(t_stats, 'spout_contact2')
                            if lick_num(k) == 1 
                                if strcmp(ylabel_input,'reaction_time') || true_L1 == 1
                                    angle_temp(m) = {prot_angle_all_session([t_stats.fakeout_trial] == j & [t_stats.lick_index] == lick_num(k) & [t_stats.time_rel_cue] > 0 & [t_stats.lesion_trial] == (m-1) & trial_num_half)};
                                else
                                    angle_temp(m) = {prot_angle_all_session([t_stats.fakeout_trial] == j & lick_index_contact == lick_num(k) & [t_stats.lesion_trial] == (m-1) & trial_num_half)};
                                end
                            else
                                angle_temp(m) = {prot_angle_all_session([t_stats.fakeout_trial] == j & ismember([t_stats.trial_num], trial_num_L2contact) & lick_index_contact == lick_num(k) & [t_stats.lesion_trial] == (m-1) & ismember([t_stats.trial_num], trial_num_L2_double_tap) & ismember([t_stats.trial_num], trial_num_L1_double_tap) & trial_num_half)};
                            end
                        else
                            disp('Warning - no spout_contact2 field in t_stats!')
                            angle_temp(m) = {prot_angle_all_session([t_stats.fakeout_trial] == j & ismember([t_stats.trial_num], trial_num_L2contact) & lick_index_contact == lick_num(k) & [t_stats.lesion_trial] == (m-1) & trial_num_half)};
                        end
                    elseif lick_num(k) == 3 % for L3 in non-recentering analyses
                        angle_temp(m) = {prot_angle_all_session([t_stats.fakeout_trial] == j & ismember([t_stats.trial_num], trial_num_L2contact) & lick_index_contact == lick_num(k) & [t_stats.lesion_trial] == (m-1) & ismember([t_stats.trial_num], trial_num_L2_double_tap) & ismember([t_stats.trial_num], trial_num_L1_double_tap) & trial_num_half )};
                    elseif lick_num(k) > 3
                        angle_temp(m) = {prot_angle_all_session([t_stats.fakeout_trial] == j & ismember([t_stats.trial_num], trial_num_L2contact) & lick_index_contact == lick_num(k) & [t_stats.lesion_trial] == (m-1) & ismember([t_stats.trial_num], trial_num_L2_double_tap) & ismember([t_stats.trial_num], trial_num_L1_double_tap) & trial_num_half)};
                    end

                    angle_temp{m} = angle_temp{m}(~isnan(angle_temp{m})); % To filter out NaNs. This pretty much only affects contact site analyses.

                    for n = 1:(numel(x_bin)-1)
                        angle_bin_temp(m, n) = sum(angle_temp{m} > x_bin(n) & angle_temp{m} <= x_bin(n+1));
                    end
                    prot_angle_bin{m, j, k} = angle_bin_temp(m, :);
                    prot_angle_temp{m, j, k} = angle_temp{m};
                end
            end
        end
    end

    % calculate median, 25th and 75th percentiles
    med_prot_angle = cellfun(@median, prot_angle_temp);
    prctile_prot_angle = cellfun(@(x) prctile((x), [25 75]), prot_angle_temp, 'UniformOutput', false);

    % delta_prot_angle organized by percentile (1 - 25, 2 - 75), condition (1 - pre-lesion, 2 - post-lesion), spout position (l/c/r, 1/2/3), lick number 
    delta_prot_angle(1,:,:,:) = med_prot_angle - cell2mat(cellfun(@(x) x(1), prctile_prot_angle, 'UniformOutput', false));
    delta_prot_angle(2,:,:,:) = cell2mat(cellfun(@(x) x(2), prctile_prot_angle, 'UniformOutput', false)) - med_prot_angle;

    % plot all in subplots
    if plot_individual == 1
        figure('Name', fig_name,'NumberTitle','off'); 
        colors = {'red', 'black', 'blue'};
        w = 1:4:20;
        x = 2:4:20;
        y = 3:4:20;
        z = 4:4:20;
        for s = 1:size(prot_angle_bin, 3)
            v = 0;
            for t = 1:size(prot_angle_bin, 2)
                % pre-lesion histogram
                subplot(numel(lick_num),4,w(s)); hold on  
                stairs(x_bin, (prot_angle_bin{1, t, s}/sum(prot_angle_bin{1, t, s})), colors{t})       
                ylim([0 0.5])
                xlim([-22 22])
                xlabel('Protrusion Angle (deg)')
                ylabel('Probability')
                title(sprintf('%s %d %s', 'Lick', s, 'Relative to First Contact (Intact)'))
                
                % post-lesion histogram
                subplot(numel(lick_num),4,x(s)); hold on  
                stairs(x_bin, (prot_angle_bin{2, t, s}/sum(prot_angle_bin{2, t, s})), colors{t})       
                ylim([0 0.5])
                xlim([-25 25])
                xlabel('Protrusion Angle (deg)')
                ylabel('Probability')
                title(sprintf('%s %d %s', 'Lick', s, 'Relative to First Contact (Inactivation)'))
                
                % median angle plot intact vs. inactivation
                subplot(numel(lick_num),4,y(s)); hold on 
                for u = 1:size(med_prot_angle, 1)
                    if u == 1
                        scatter(u+v, med_prot_angle(u, t, s), [], colors{t}) 
                        errorbar(u+v,  med_prot_angle(u, t, s), delta_prot_angle(1,u,t,s), delta_prot_angle(2,u,t,s), 'color', colors{t})
                    elseif u == 2
                        scatter(u+v, med_prot_angle(u, t, s), [], [0.3010 0.7450 0.9330]) 
                        errorbar(u+v,  med_prot_angle(u, t, s), delta_prot_angle(1,u,t,s), delta_prot_angle(2,u,t,s), 'color', [0.3010 0.7450 0.9330])
                    end
                end
                xlim([0 9])
                ylim(ylim_val)
                xticks(1.5:3:7.5)
                xticklabels({'Left', 'Center', 'Right'})
                xlabel('Spout Position')
                ylabel('Angle (deg)')
                title(sprintf('%s %d %s', 'Lick', s, 'Relative to First Contact'))                
                v = v + 3; % increase plotting index for x-axis
            end
        end
    end

    med_prot_angle_all(p,:,:,:) = num2cell(med_prot_angle);
    prot_angle(p, :, :, :) = prot_angle_temp;

    % concatenate across conditions for permutation test
    for q = 1:size(prot_angle, 3)
        for r = 1:size(prot_angle, 4)
            prot_angle_all_lesion{p, q, r} = vertcat(prot_angle_temp{:, q, r});
        end
    end
    
    % concatenate across positions for permutation test
    for q = 1:size(prot_angle, 2)
        for r = 1:size(prot_angle, 4)
            prot_angle_all_position{p, q, r} = vertcat(prot_angle_temp{q, :, r});
        end
    end    
end

% clearvars -except prot_angle prot_angle_all_laser prot_angle_all_position prot_lat_disp prot_lat_disp_all_laser prot_lat_disp_all_position ylim_val ylabel_input difference_laser correction_sum opto_side correction_sum_xsides correction_mag_sum all_values prctile_prot_angle

% calculate true medians and difference across population of mice
med_prot_angle_temp = cellfun(@median, prot_angle);
med_prot_angle_all = squeeze(median(med_prot_angle_temp, 1));
med_prot_angle_diff = abs(squeeze(med_prot_angle_all(1,:,:)) - squeeze(med_prot_angle_all(2,:,:)));

numShuffles = 1000; % perform a hierarchical bootstrap to get percentiles

% shuf_prot_angle organized by number of shuffles, animal, condition (1 - pre-lesion, 2 - post-lesion), spout position (1/2/3, left/center/right), lick number
shuf_med_prot_angle = zeros([numShuffles size(prot_angle)]);

for i = 1:numShuffles

    % randomly select mice to create a new pseudopopulation
    shuf_animal_ind = randi(size(prot_angle, 1), 1, size(prot_angle, 1));
    shuf_angle_animal = prot_angle(shuf_animal_ind, :, :, :);
    
    % randomly select trials from each mouse/condition to get a new pseudosession
    for j = 1:size(shuf_angle_animal, 1)
        for k = 1:size(shuf_angle_animal, 2)
            trial_num_temp = cellfun(@numel, squeeze(shuf_angle_animal(j, k, :, :)), 'UniformOutput', false);
            
            if sum([trial_num_temp{:}]) == 0
                shuf_med_prot_angle(i, j, k, :, :) = [NaN;NaN;NaN];
            else
                trial_num_ind_temp = cellfun(@(x) randi(x, 1, x), trial_num_temp, 'UniformOutput', false);
                shuf_prot_angle_temp = cellfun(@(y, z) y(z), squeeze(shuf_angle_animal(j, k, :, :)), trial_num_ind_temp, 'UniformOutput', false);
                shuf_med_prot_angle(i, j, k, :, :) = cellfun(@median, shuf_prot_angle_temp);
            end
        end
    end
end

% calculate median and interquartile range
shuf_med_prot_angle_all = squeeze(median(shuf_med_prot_angle, 2));
shuf_med_prot_angle_all_prctile = squeeze(prctile(shuf_med_prot_angle_all, [25 75], 1));

% calculate deltas for errorbar plotting - array organized by errobar (1 - lower, 2 - upper), condition (1 - pre-lesion, 2 - post-lesion), spout position(1/2/3, left/center/right), lick number
delta_prot_angle_all(1,:,:,:) = med_prot_angle_all - squeeze(shuf_med_prot_angle_all_prctile(1,:,:,:));
delta_prot_angle_all(2,:,:,:) = squeeze(shuf_med_prot_angle_all_prctile(2,:,:,:)) - med_prot_angle_all;

shuffle_test_lesion_tailOne = nan(numShuffles,size(shuf_med_prot_angle_all,3));
shuffle_test_lesion_tailTwo = nan(numShuffles,size(shuf_med_prot_angle_all,3));
shuffle_test_position_tailOne = nan(numShuffles,size(shuf_med_prot_angle_all,2),size(shuf_med_prot_angle_all,3));
shuffle_test_position_tailTwo = nan(numShuffles,size(shuf_med_prot_angle_all,2),size(shuf_med_prot_angle_all,3));
for i = 1:numShuffles
    shuffle_test_lesion_tailOne(i,:) = shuf_med_prot_angle_all(i,1,:) - shuf_med_prot_angle_all(i,2,:) <= 0;
    shuffle_test_lesion_tailTwo(i,:) = shuf_med_prot_angle_all(i,1,:) - shuf_med_prot_angle_all(i,2,:) >= 0;
    for m = 1:size(shuf_med_prot_angle_all,2)
        for n = 1:size(shuf_med_prot_angle_all,3)
            shuffle_test_position_tailOne(i,m,n) = shuf_med_prot_angle_all(i,m,n) - shuf_med_prot_angle_all(i,m,2) <= 0;
            shuffle_test_position_tailTwo(i,m,n) = shuf_med_prot_angle_all(i,m,n) - shuf_med_prot_angle_all(i,m,2) >= 0;
        end      
    end
end
shuffle_test_lesion_tailOne = mean(shuffle_test_lesion_tailOne,1);
shuffle_test_lesion_tailTwo = mean(shuffle_test_lesion_tailTwo,1);
shuffle_test_position_tailOne = squeeze(mean(shuffle_test_position_tailOne,1));
shuffle_test_position_tailTwo = squeeze(mean(shuffle_test_position_tailTwo,1));

% plot medians + bootstrapped IQR for across all mice
plot_data_indiv = [];
plot_data_summary = [];
figure('Name', fig_name,'NumberTitle','off'); 
w = 1:2:6;
med_prot_angle = cellfun(@median, prot_angle);
for s = 1:size(med_prot_angle_all, 3)
    v = 0;
    for t = 1:size(med_prot_angle_all, 2)
        % median angle plot intact vs. inactivation
        subplot(size(med_prot_angle_all, 3),2,w(s)); hold on 
        for u = 1:size(med_prot_angle_all, 1)
            if u == 1
                scatter(u+v, med_prot_angle_all(u, t, s), [], [0 0 0],"filled") 
                errorbar(u+v,  med_prot_angle_all(u, t, s), delta_prot_angle_all(1,u,t,s), delta_prot_angle_all(2,u,t,s), 'LineWidth', 1, 'color', [0.5 0.5 0.5])
            elseif u == 2
                scatter(u+v, med_prot_angle_all(u, t, s), [], [0.3010 0.7450 0.9330],"filled") 
                errorbar(u+v,  med_prot_angle_all(u, t, s), delta_prot_angle_all(1,u,t,s), delta_prot_angle_all(2,u,t,s), 'LineWidth', 1, 'color', [0.5 0.5 0.5])
            end
            plot_data_summary = [plot_data_summary;[u+v,  med_prot_angle_all(u, t, s), med_prot_angle_all(u, t, s)-delta_prot_angle_all(1,u,t,s), med_prot_angle_all(u, t, s)+delta_prot_angle_all(2,u,t,s)]];
            temp_array=med_prot_angle(:,:,t);
            scatter(u+v, temp_array(:,u), [], [.7 .7 .7],"filled")
            plot_data_indiv = [plot_data_indiv;[u+v, [temp_array(:,u)]']];
        end
        xlim([0 9])
        ylim(ylim_val)
        xticks(1.5:3:7.5)
        xticklabels({'Left', 'Center', 'Right'})
        xlabel('Spout Position')
        ylabel(ylabel_input)
        title(sprintf('%s %d %s', 'Lick', s, 'Relative to First Contact'))
        v = v + 3;
    end
end

% makes it easier to report IQR
iqr_prot_anlgle_all = cell(size(med_prot_angle_all,1),size(med_prot_angle_all,2));
for i = 1:size(med_prot_angle_all,1)
    for j = 1: size(med_prot_angle_all,2)
        iqr_prot_anlgle_all{i,j}{1} = med_prot_angle_all(i,j)-delta_prot_angle_all(1,i,j);
        iqr_prot_anlgle_all{i,j}{2} = med_prot_angle_all(i,j)+delta_prot_angle_all(2,i,j);
    end
end

% the following code combines data across trial types
prot_angle_combined = cell(size(prot_angle,1),size(prot_angle,2));
for i = 1:size(prot_angle,1)
    for j = 1:size(prot_angle,2)
        prot_angle_combined{i,j} = [prot_angle{i,j,1};prot_angle{i,j,2};prot_angle{i,j,3}];
    end
end

med_prot_angle_temp_combined = cellfun(@median, prot_angle_combined);
med_prot_angle_all_combined = squeeze(median(med_prot_angle_temp_combined, 1));
shuf_med_prot_angle_combined = zeros([numShuffles size(prot_angle_combined)]);
for i = 1:numShuffles

    % randomly select mice to create a new pseudopopulation
    shuf_animal_ind = randi(size(prot_angle_combined, 1), 1, size(prot_angle_combined, 1));
    shuf_angle_animal = prot_angle_combined(shuf_animal_ind, :);
    
    % randomly select trials from each mouse/condition to get a new pseudosession
    for j = 1:size(shuf_angle_animal, 1)
        for k = 1:size(shuf_angle_animal, 2)
            trial_num_temp = cellfun(@numel, squeeze(shuf_angle_animal(j, k, :, :)), 'UniformOutput', false);
            
            if sum([trial_num_temp{:}]) == 0
                shuf_med_prot_angle_combined(i, j, k) = [NaN;NaN;NaN];
            else
                trial_num_ind_temp = cellfun(@(x) randi(x, 1, x), trial_num_temp, 'UniformOutput', false);
                shuf_prot_angle_temp = cellfun(@(y, z) y(z), squeeze(shuf_angle_animal(j, k, :, :)), trial_num_ind_temp, 'UniformOutput', false);
                shuf_med_prot_angle_combined(i, j, k) = cellfun(@median, shuf_prot_angle_temp);
            end
        end
    end
end

shuf_med_prot_angle_combined_all = squeeze(median(shuf_med_prot_angle_combined, 2));

shuffle_test_lesion_combined_tailOne = nan(numShuffles,1);
shuffle_test_lesion_combined_tailTwo = nan(numShuffles,1);
for i = 1:numShuffles
    shuffle_test_lesion_combined_tailOne(i) = shuf_med_prot_angle_combined_all(i,1) - shuf_med_prot_angle_combined_all(i,2) <= 0;
    shuffle_test_lesion_combined_tailTwo(i) = shuf_med_prot_angle_combined_all(i,1) - shuf_med_prot_angle_combined_all(i,2) >= 0;
end
shuffle_test_lesion_combined_tailOne = mean(shuffle_test_lesion_combined_tailOne,1);
shuffle_test_lesion_combined_tailTwo = mean(shuffle_test_lesion_combined_tailTwo,1);

figure('Name', fig_name,'NumberTitle','off'); 
hold on
med_prot_angle_combined = cellfun(@median, prot_angle_combined);
for u = 1:size(med_prot_angle_all_combined,2)
    if u == 1
        scatter(1.2*ones(1,size(med_prot_angle_combined,1)), med_prot_angle_combined(:,u), [], [0 0 0],"filled")
        errorbar(1.2,  med_prot_angle_all_combined(u), med_prot_angle_all_combined(u)-prctile(shuf_med_prot_angle_combined_all(:,u),25), prctile(shuf_med_prot_angle_combined_all(:,u),75)-med_prot_angle_all_combined(u), 'LineWidth', 10, 'color', [0.5 0.5 0.5])
    elseif u == 2
        scatter(1.8*ones(1,size(med_prot_angle_combined,1)), med_prot_angle_combined(:,u), [], [0.3010 0.7450 0.9330],"filled") 
        errorbar(1.8,  med_prot_angle_all_combined(u), med_prot_angle_all_combined(u)-prctile(shuf_med_prot_angle_combined_all(:,u),25), prctile(shuf_med_prot_angle_combined_all(:,u),75)-med_prot_angle_all_combined(u), 'LineWidth', 10, 'color', [0.5 0.5 0.5])
    end
end

plot_data_summary_combined = [med_prot_angle_all_combined(1), prctile(shuf_med_prot_angle_combined_all(:,1),25), prctile(shuf_med_prot_angle_combined_all(:,1),75);...
                     med_prot_angle_all_combined(2), prctile(shuf_med_prot_angle_combined_all(:,2),25), prctile(shuf_med_prot_angle_combined_all(:,2),75)];
plot_data_indiv_combined = [med_prot_angle_combined(:,1);med_prot_angle_combined(:,2)];

xlim([1 2])
ylim(ylim_val)
xticks([1.2,1.8])
xticklabels({'Pre-lesion', 'Post-lesion'})
ylabel(ylabel_input)
if lick_num == 1
    title(sprintf('%s %d %s', 'First Lick'))
else
    title(sprintf('%s %d %s', ['Lick ', num2str(lick_num), ' Relative to First Contact']))
end
IQR_for_table = [prctile(shuf_med_prot_angle_combined_all(:,1),25), prctile(shuf_med_prot_angle_combined_all(:,1),75);prctile(shuf_med_prot_angle_combined_all(:,2),25), prctile(shuf_med_prot_angle_combined_all(:,2),75)];

end