function [plot_data_summary,plot_data_indiv,shuffle_test_laser_tailOne,shuffle_test_laser_tailTwo] = plot_probability_trial_type(data_set,trial_type,probability_type, fig_name)

for k = 1:numel(data_set)
    load(data_set{k});
    second_half = k-1; % if this iteration is for the second unilateral side, combine data with the previous side 
    reverse_index = [numel(trial_type):-1:1];

    if second_half == 0
        prot_num_first_half = zeros(numel(t_stats_combined),2,numel(trial_type));
        trial_num_first_half = zeros(numel(t_stats_combined),2,numel(trial_type));
        raw_data = cell (numel(t_stats_combined),2,numel(trial_type));
        prot_prob = [];
    end

    for i = 1:numel(t_stats_combined)
        t_stats = t_stats_combined(i).data;   
        for j = 1:numel(trial_type)
            trial_num_L1_double_tap = [t_stats(isnan([t_stats.spout_contact2]) & [t_stats.lick_index_contact1] == 1).trial_num];
            exclude_trials = [t_stats(~isnan([t_stats.spout_contact2]) & isnan([t_stats.spout_contact])).trial_num]; % sometimes a lick would have no spoutcontact but do have spoutcontact2, and those are trials that mess up with the lick rel2 contact indices

            if probability_type == 1
                %for contact probability, only calculate probablility for trials where there is an L3
                index_lick1 = find(ismember([t_stats.fakeout_trial],trial_type(j)) & ismember([t_stats.trial_num],trial_num_L1_double_tap) & ~ismember([t_stats.trial_num],exclude_trials) & ([t_stats.lick_index_contact2] == 2) & ([t_stats.laser_trial] == 0) & ~isnan([t_stats.spout_contact]) & ~isnan([t_stats.prev_spcontact]));
                index_lick2 = find(ismember([t_stats.fakeout_trial],trial_type(j)) & ismember([t_stats.trial_num],trial_num_L1_double_tap) & ~ismember([t_stats.trial_num],exclude_trials) & ([t_stats.lick_index_contact2] == 2) & ([t_stats.laser_2D] == 1) & ~isnan([t_stats.spout_contact]) & ~isnan([t_stats.prev_spcontact]));
                prot_nl = [t_stats(index_lick1).prev_spcontact];
                prot_l = [t_stats(index_lick2).prev_spcontact];
                num_control_all = numel(find(ismember([t_stats.fakeout_trial],trial_type(j)) & ismember([t_stats.trial_num],trial_num_L1_double_tap) & ~ismember([t_stats.trial_num],exclude_trials) ...
                               & [t_stats.lick_index_contact2] == 2 & [t_stats.laser_trial] == 0 & ~isnan([t_stats.prev_spcontact]))); % don't count trials where L1 already triggered laser
                num_laser_all = numel(find(ismember([t_stats.fakeout_trial],trial_type(j)) & ismember([t_stats.trial_num],trial_num_L1_double_tap) & ~ismember([t_stats.trial_num],exclude_trials)...
                               & [t_stats.lick_index_contact2] == 2 & [t_stats.laser_2D] == 1 & ~isnan([t_stats.prev_spcontact])));
                if second_half == 1 % if it's the second half, append to the existing data. Additionally, flip left and right trials so that ipsi/contra sides match between the two unilateral inactivations.
                    if i <= size(prot_num_first_half,1)
                        prot_prob(i,1,reverse_index(j)) = (numel(find(prot_nl<=250)) + prot_num_first_half(i,1,reverse_index(j))) / (num_control_all + trial_num_first_half(i,1,reverse_index(j)));
                        prot_prob(i,2,reverse_index(j)) = (numel(find(prot_l<=250)) + prot_num_first_half(i,2,reverse_index(j)))  / (num_laser_all + trial_num_first_half(i,2,reverse_index(j)));
                    else
                        prot_prob(i,1,reverse_index(j)) = numel(find(prot_nl<=250)) / num_control_all;
                        prot_prob(i,2,reverse_index(j)) = numel(find(prot_l<=250)) / num_laser_all;
                    end
                elseif second_half == 0
                        prot_prob(i,1,j) = numel(find(prot_nl<=250)) / num_control_all;
                        prot_prob(i,2,j) = numel(find(prot_l<=250)) / num_laser_all;
                end
            elseif probability_type == 2
                % for CSM probability, only calculate probablility for trials where there is an L3
                index_lick1 = find(ismember([t_stats.fakeout_trial],trial_type(j)) & ismember([t_stats.trial_num],trial_num_L1_double_tap) & ~ismember([t_stats.trial_num],exclude_trials) & ([t_stats.lick_index_contact2] == 2) & ([t_stats.laser_trial] == 0));
                index_lick2 = find(ismember([t_stats.fakeout_trial],trial_type(j)) & ismember([t_stats.trial_num],trial_num_L1_double_tap) & ~ismember([t_stats.trial_num],exclude_trials) & ([t_stats.lick_index_contact2] == 2) & ([t_stats.laser_2D] == 1) );
                prot_nl = [t_stats(index_lick1).CSM_dur];
                prot_l = [t_stats(index_lick2).CSM_dur];
                num_control_all = numel(find(ismember([t_stats.fakeout_trial],trial_type(j)) & ismember([t_stats.trial_num],trial_num_L1_double_tap) & ~ismember([t_stats.trial_num],exclude_trials) ...
                                  & [t_stats.lick_index_contact2] == 2 & [t_stats.laser_trial] == 0));        
                num_laser_all = numel(find(ismember([t_stats.fakeout_trial],trial_type(j)) & ismember([t_stats.trial_num],trial_num_L1_double_tap) & ~ismember([t_stats.trial_num],exclude_trials) & [t_stats.lick_index_contact2] == 2 & [t_stats.laser_trial] == 1)); % don't count trials where L1 already triggered laser
                if second_half == 1 % if it's the second half, append to the existing data. Additionally, flip left and right trials so that ipsi/contra sides match between the two unilateral inactivations.
                    if i <= size(prot_num_first_half,1)
                        prot_prob(i,1,reverse_index(j)) = (numel(find(prot_nl>5)) + prot_num_first_half(i,1,reverse_index(j))) / (num_control_all + trial_num_first_half(i,1,reverse_index(j)));
                        prot_prob(i,2,reverse_index(j)) = (numel(find(prot_l>5)) + prot_num_first_half(i,2,reverse_index(j)))  / (num_laser_all + trial_num_first_half(i,2,reverse_index(j)));
                    else
                        prot_prob(i,1,reverse_index(j)) = numel(find(prot_nl>5)) / num_control_all;
                        prot_prob(i,2,reverse_index(j)) = numel(find(prot_l>5)) / num_laser_all;
                    end
                elseif second_half == 0
                    prot_prob(i,1,j) = numel(find(prot_nl>5)) / num_control_all;
                    prot_prob(i,2,j) = numel(find(prot_l>5)) / num_laser_all;
                end
            end
            
            if second_half == 0        
                if probability_type == 1
                    prot_num_first_half(i,1,j) = numel(find(prot_nl<=250));
                    prot_num_first_half(i,2,j) = numel(find(prot_l<=250));                
                elseif probability_type == 2
                    prot_num_first_half(i,1,j) = numel(find(prot_nl>5));
                    prot_num_first_half(i,2,j) = numel(find(prot_l>5));
                end
                trial_num_first_half(i,1,j) = num_control_all;
                trial_num_first_half(i,2,j) = num_laser_all;
                raw_data{i,1,j} = vertcat(prot_nl',nan(num_control_all-numel(prot_nl),1));
                raw_data{i,2,j} = vertcat(prot_l',nan(num_laser_all-numel(prot_l),1));
            elseif second_half == 1 && i <= size(raw_data,1)
                raw_data{i,1,reverse_index(j)} = vertcat(raw_data{i,1,reverse_index(j)}, prot_nl',nan(num_control_all-numel(prot_nl),1));
                raw_data{i,2,reverse_index(j)} = vertcat(raw_data{i,2,reverse_index(j)}, prot_l',nan(num_laser_all-numel(prot_l),1));
            else
                raw_data{i,1,reverse_index(j)} = vertcat(prot_nl',nan(num_control_all-numel(prot_nl),1));
                raw_data{i,2,reverse_index(j)} = vertcat(prot_l',nan(num_laser_all-numel(prot_l),1));
            end
        end
    end
end

figure('Name', fig_name,'NumberTitle','off');
colors = {'red', 'black', 'blue'};
w = 1:2:6;
x = 2:2:6;

loop_num = 1000;
perm_med_prot_angle_laser_tailOne = nan([loop_num size(prot_prob,3)]);
perm_med_prot_angle_laser_tailTwo = nan([loop_num size(prot_prob,3)]);
all_medians=nan(size(prot_prob,2),size(prot_prob,3),loop_num);
med_prob_all = nan(size(prot_prob,2),size(prot_prob,3));
delta_prob_all = nan(2,size(prot_prob,2),size(prot_prob,3));

for k=1:loop_num
    ID_BS_animal = datasample([1:size(raw_data,1)],size(raw_data,1));
    data_temp = raw_data(ID_BS_animal,:,:);
    data_BS_temp = cell(size(data_temp,2),size(data_temp,3));
    
    for i = 1:size(data_temp,3)
        for j = 1:size(data_temp,2)
            data_BS_temp_condition = data_temp(:,j,i);
            probability_BS_temp=nan(1,numel(data_BS_temp_condition));
            for l=1:numel(data_BS_temp_condition)
                temp=datasample(data_BS_temp_condition{l},numel(data_BS_temp_condition{l}));                
                if probability_type == 1
                    probability_BS_temp(l)=numel(find(temp<=250))/numel(data_BS_temp_condition{l});
                elseif probability_type == 2                
                    probability_BS_temp(l)=numel(find(temp>5))/numel(data_BS_temp_condition{l});
                end
            end
            data_BS_temp{j,i} = probability_BS_temp;
            all_medians(j,i,k) = median(probability_BS_temp);
        end
        perm_med_prot_angle_laser_tailOne(k,i) = median(data_BS_temp{1,i}) - median(data_BS_temp{2,i}) <= 0;
        perm_med_prot_angle_laser_tailTwo(k,i) = median(data_BS_temp{1,i}) - median(data_BS_temp{2,i}) >= 0;
    end
end

shuffle_test_laser_tailOne = mean(perm_med_prot_angle_laser_tailOne,1);
shuffle_test_laser_tailTwo = mean(perm_med_prot_angle_laser_tailTwo,1);

for i = 1:size(prot_prob,2)
    for j = 1:size(prot_prob,3) 
        med_prob_all(i,j) = median(prot_prob(:,i,j));
        delta_prob_all(1,i,j) = prctile(all_medians(i,j,:),25);
        delta_prob_all(2,i,j) = prctile(all_medians(i,j,:),75);        
    end
end

plot_data_indiv = [];
plot_data_summary = [];
for s = 1:size(med_prob_all, 3)
    v = 0;
    for t = 1:size(med_prob_all, 2)
        % median angle plot intact vs. inactivation
        subplot(size(med_prob_all, 3),2,w(s)); hold on 
        for u = 1:size(med_prob_all, 1)
            if u == 1
                scatter(u+v, med_prob_all(u, t, s), [], [0 0 0],"filled") 
                errorbar(u+v,  med_prob_all(u, t, s), med_prob_all(u, t, s)-delta_prob_all(1,u,t,s), delta_prob_all(2,u,t,s)-med_prob_all(u, t, s), 'LineWidth', 1, 'color', [0 0 0])
            elseif u == 2
                scatter(u+v, med_prob_all(u, t, s), [], [0.3010 0.7450 0.9330],"filled")
                errorbar(u+v,  med_prob_all(u, t, s), med_prob_all(u, t, s)-delta_prob_all(1,u,t,s), delta_prob_all(2,u,t,s)-med_prob_all(u, t, s), 'LineWidth', 1, 'color', [0.3010 0.7450 0.9330])
            end
            plot_data_summary = [plot_data_summary;[u+v,  med_prob_all(u, t, s), med_prob_all(u, t, s)-delta_prob_all(1,u,t,s), med_prob_all(u, t, s)+delta_prob_all(2,u,t,s)]];
            temp_array=prot_prob(:,:,t);
            scatter(u+v, temp_array(:,u), [], [.7 .7 .7],"filled")
            plot_data_indiv = [plot_data_indiv;[u+v, [temp_array(:,u)]']];
        end
        xlim([0 9])
        ylim([-0.05 1.05])
        xticks(1.5:3:7.5)
        xticklabels({'Left', 'Center', 'Right'})
        xlabel('Spout Position')
        ylabel('L3 Prot. Probability')
        title(sprintf('%s %d %s', 'Lick', s, 'Relative to First Contact'))
        v = v + 3;
    end
end
end