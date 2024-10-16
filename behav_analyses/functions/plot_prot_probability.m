function [plot_data_summary,plot_data_indiv,shuffle_test_laser_tailOne,shuffle_test_laser_tailTwo] = plot_prot_probability(data_set,trial_type,probability_type, fig_name)

load(data_set)
num_laser_all = 0;
trial_num_all = 0;
prot_nl_all = [];
prot_l_all = [];
prot_nl_byMouse = cell(1,numel(t_stats_combined));
prot_l_byMouse = cell(1,numel(t_stats_combined));

for i = 1:numel(t_stats_combined)
    t_stats = t_stats_combined(i).data;   
    trial_num_L1_double_tap = [t_stats(isnan([t_stats.spout_contact2]) & [t_stats.lick_index_contact1] == 1).trial_num];
    exclude_trials = [t_stats(~isnan([t_stats.spout_contact2]) & isnan([t_stats.spout_contact])).trial_num]; % sometimes a lick would have no spoutcontact but do have spoutcontact2, and those are trials that mess up with the lick rel2 contact indices

    if probability_type == 1
        index_lick1 = find(ismember([t_stats.fakeout_trial],trial_type) & ismember([t_stats.trial_num],trial_num_L1_double_tap) & ~ismember([t_stats.trial_num],exclude_trials) & ([t_stats.lick_index_contact2] == 2) & ([t_stats.laser_trial] == 0)& ~isnan([t_stats.prev_spcontact]));
        index_lick2 = find(ismember([t_stats.fakeout_trial],trial_type) & ismember([t_stats.trial_num],trial_num_L1_double_tap) & ~ismember([t_stats.trial_num],exclude_trials) & ([t_stats.lick_index_contact2] == 2) & ([t_stats.laser_2D] == 1)& ~isnan([t_stats.prev_spcontact]));
        prot_nl = [t_stats(index_lick1).time_rel_cue] - [t_stats(index_lick1-1).spout_contact];
        prot_l = [t_stats(index_lick2).time_rel_cue] - [t_stats(index_lick2-1).spout_contact];
    elseif probability_type == 2 % for contact probability, only considers L3 that made contacts
        index_lick1 = find(ismember([t_stats.fakeout_trial],trial_type) & ismember([t_stats.trial_num],trial_num_L1_double_tap) & ~ismember([t_stats.trial_num],exclude_trials) & ([t_stats.lick_index_contact2] == 2) & ([t_stats.laser_trial] == 0) & ~isnan([t_stats.prev_spcontact]) & ~isnan([t_stats.spout_contact]));
        index_lick2 = find(ismember([t_stats.fakeout_trial],trial_type) & ismember([t_stats.trial_num],trial_num_L1_double_tap) & ~ismember([t_stats.trial_num],exclude_trials) & ([t_stats.lick_index_contact2] == 2) & ([t_stats.laser_2D] == 1) & ~isnan([t_stats.prev_spcontact]) & ~isnan([t_stats.spout_contact]));
        prot_nl = [t_stats(index_lick1).prev_spcontact];
        prot_l = [t_stats(index_lick2).prev_spcontact];
    end

    % count trial numbers
    laser_num_temp = numel(find(ismember([t_stats.fakeout_trial],trial_type) & ismember([t_stats.trial_num],trial_num_L1_double_tap) & ~ismember([t_stats.trial_num],exclude_trials) & [t_stats.lick_index_contact2] == 1 & [t_stats.laser_trial] == 1)); % don't count trials where L1 already triggered laser
    num_laser_all = num_laser_all + laser_num_temp; 
    all_num_temp = numel(find(ismember([t_stats.fakeout_trial],trial_type) & ismember([t_stats.trial_num],trial_num_L1_double_tap) & ~ismember([t_stats.trial_num],exclude_trials) & [t_stats.lick_index_contact2] == 1 ));
    trial_num_all = trial_num_all + all_num_temp;
    % collect all latencies
    prot_nl_all = [prot_nl_all prot_nl];
    prot_l_all = [prot_l_all prot_l];
    % collect latencies for each mouse
    prot_nl_byMouse{i} = [prot_nl, nan(1,all_num_temp - laser_num_temp - numel(prot_nl))];
    prot_l_byMouse{i} = [prot_l, nan(1,laser_num_temp - numel(prot_l))];
end

median_prot_nl(i) = median(prot_nl_all);
median_prot_l(i) = median(prot_l_all);

% plot culmulative probability curves
edges = 0:1:1000;
prot_nl_vect = histc(prot_nl_all,edges);
prot_l_vect = histc(prot_l_all,edges);
plot_vect_nl(i,:) = cumsum(prot_nl_vect)/(trial_num_all-num_laser_all);

if numel(prot_l)
    plot_vect_l(i,:) = cumsum(prot_l_vect)/num_laser_all;
else
    plot_vect_l(i,:) = zeros(1,numel(edges));
end

figure('Name', fig_name,'NumberTitle','off'); 
plot(plot_vect_nl(i,:),'k');hold on;
plot(plot_vect_l(i,:),'b');
xlim([1 250]);
ylim([-0.1 1]);
set(gcf,'position',[560 210 400 400])

% get and plot probability of laser vs no laser conditions
loop_num = 1000;
median_nl = median(cellfun(@(x) sum(x<=250)/numel(x), prot_nl_byMouse));
medial_l = median(cellfun(@(x) sum(x<=250)/numel(x), prot_l_byMouse));
median_nl_BS = nan(1,loop_num);
median_l_BS = nan(1,loop_num);
shuffle_test_tailOne = nan(1,loop_num);
shuffle_test_tailTwo = nan(1,loop_num);

for i = 1:loop_num
    bs_animal = datasample([1:numel(t_stats_combined)],numel(t_stats_combined));
    nl_temp = prot_nl_byMouse(bs_animal);
    l_temp = prot_l_byMouse(bs_animal);
    for j = 1:numel(nl_temp)
        nl_temp{j} = nl_temp{j}(datasample([1:numel(nl_temp{j})], numel(nl_temp{j})));
        l_temp{j} = l_temp{j}(datasample([1:numel(l_temp{j})], numel(l_temp{j})));
    end
    median_nl_BS(i) = median(cellfun(@(x) sum(x<=250)/numel(x), nl_temp));    
    median_l_BS(i) = median(cellfun(@(x) sum(x<=250)/numel(x), l_temp));
    shuffle_test_tailOne(i) = median_nl_BS(i) - median_l_BS(i) <=0;
    shuffle_test_tailTwo(i) = median_nl_BS(i) - median_l_BS(i) >=0;
end

shuffle_test_laser_tailOne = mean(shuffle_test_tailOne);
shuffle_test_laser_tailTwo = mean(shuffle_test_tailTwo);

figure('Name', fig_name,'NumberTitle','off'); 
hold on
scatter(1.2*ones(1,numel(prot_nl_byMouse)), cellfun(@(x) sum(x<=250)/numel(x), prot_nl_byMouse), [], [0.7 0.7 0.7],"filled")
errorbar(1.2,median_nl, median_nl-prctile(median_nl_BS,25), prctile(median_nl_BS,75)-median_nl, 'LineWidth', 1, 'color', [0.3010 0.7450 0.9330],'Marker','.','MarkerSize',12)
scatter(1.8*ones(1,numel(prot_l_byMouse)), cellfun(@(x) sum(x<=250)/numel(x), prot_l_byMouse), [], [0.7 0.7 0.7],"filled")
errorbar(1.8,medial_l, medial_l-prctile(median_l_BS,25), prctile(median_l_BS,75)-medial_l, 'LineWidth', 1, 'color', [0.3010 0.7450 0.9330],'Marker','.','MarkerSize',12)
plot_data_summary = [median_nl, prctile(median_nl_BS,25), prctile(median_nl_BS,75);...
                     medial_l, prctile(median_l_BS,25), prctile(median_l_BS,75)];
plot_data_indiv = [cellfun(@(x) sum(x<=250)/numel(x), prot_nl_byMouse);...
                   cellfun(@(x) sum(x<=250)/numel(x), prot_l_byMouse)];
xlim([0.8 2.2])
ylim([-0.05 1.05])
end