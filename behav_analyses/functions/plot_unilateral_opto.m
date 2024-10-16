function [plot_data_summary,plot_data_indiv,shuffle_test_laser_tailOne,shuffle_test_laser_tailTwo] = plot_unilateral_opto(angle_sum,ylabel_input,ylim_val, fig_name)
% set the following parameter
reverse_sign = -1; % -1 (for angles) or 1 (non signed variables)
% this code assumes the first half of ang_sum is unilateral inactivations on one hemisphere, the second half on the other hemisphere
summary_data = angle_sum;

% organize data as ipsi/center/right, instead of left/center/reight, reverse +/- sign if the variable is an angle and the inactivation is on the right hemisphere
reverse_index=[size(summary_data,2):-1:1];
for i = (size(summary_data,1)/2) + 1 : size(summary_data,1)
    data_temp = {summary_data{i,:}};
    for j = 1:size(summary_data,2)
        summary_data{i,reverse_index(j)} = cellfun(@(x) x*(reverse_sign), data_temp{1,j},'un',0);
    end
end

% Combine data from the same mouse, note that some mice may only have unilateral inactivation on one hemisphere. The following code takes care of that
summary_data_median = [];
for i = 1:size(summary_data,2)
    for j = 1:size(summary_data{1,i},2)
        for k = 1:min(size(summary_data{1,i},1),size(summary_data{2,i},1))
            summary_data{1,i}{k,j} = vertcat(summary_data{1,i}{k,j}, summary_data{2,i}{k,j});
        end
        if size(summary_data{1,i},1) < size(summary_data{2,i},1)
            summary_data{1,i} = vertcat(summary_data{1,i},summary_data{2,i}(k+1:end,:));
        end
    end
    summary_data{2,i}=cellfun(@median,summary_data{1,i});
    summary_data_median(1,i) = median(summary_data{2,i}(:,1));
    summary_data_median(2,i) = median(summary_data{2,i}(:,2));
end

% perform boostrapping for IQR and stats test
loop_num = 1000;
perm_med_prot_angle_laser_tailOne = nan([loop_num size(angle_sum,2)]);
perm_med_prot_angle_laser_tailTwo = nan([loop_num size(angle_sum,2)]);
all_medians=nan(size(angle_sum,1),size(angle_sum,2),loop_num);
summary_data_delta = nan(2,size(angle_sum,1),size(angle_sum,2));

% rearange the data to make bootstrapping easier
rearranged_data = {};
for l = 1:size(summary_data{1,1},1)
    temp = vertcat(summary_data{1,1}(l,:),summary_data{1,2}(l,:),summary_data{1,3}(l,:))';
    rearranged_data = cat(3,rearranged_data,temp);
end
rearranged_data = permute(rearranged_data,[3,1,2]);

for k=1:loop_num
    ID_BS_animal = datasample([1:size(rearranged_data,1)],size(rearranged_data,1)); % resample animals with replacement
    data_temp = rearranged_data(ID_BS_animal,:,:);
    data_BS_temp = cell(size(data_temp,2),size(data_temp,3));
    
    for i = 1:size(data_temp,3)
        for j = 1:size(data_temp,2)
            data_BS_temp_condition = data_temp(:,j,i);
            median_BS_temp = nan(1,numel(data_BS_temp_condition));
            for l = 1:numel(data_BS_temp_condition)
                temp = datasample(data_BS_temp_condition{l},numel(data_BS_temp_condition{l})); % resample trials within each condition with replacement
                median_BS_temp(l) = median(temp);
            end
            data_BS_temp{j,i} = median_BS_temp;
            all_medians(j,i,k) = median(median_BS_temp);
        end
        % for each iteration, compare the median value in the laser off condition with that in the laser on condition, on both tails
        perm_med_prot_angle_laser_tailOne(k,i) = median(data_BS_temp{1,i}) - median(data_BS_temp{2,i}) <= 0;
        perm_med_prot_angle_laser_tailTwo(k,i) = median(data_BS_temp{1,i}) - median(data_BS_temp{2,i}) >= 0;
    end
end

% calculate raw p-values for laser on vs off comparisons on both tails
shuffle_test_laser_tailOne = nan(1,size(perm_med_prot_angle_laser_tailOne,2));
shuffle_test_laser_tailTwo = nan(1,size(perm_med_prot_angle_laser_tailTwo,2));
for i = 1:numel(shuffle_test_laser_tailOne)
    shuffle_test_laser_tailOne(i) = sum(perm_med_prot_angle_laser_tailOne(:,i))/size(perm_med_prot_angle_laser_tailOne,1);
    shuffle_test_laser_tailTwo(i) = sum(perm_med_prot_angle_laser_tailTwo(:,i))/size(perm_med_prot_angle_laser_tailTwo,1);
end

% calculate 25th and 75th percentiles
for i = 1:size(angle_sum,1)
    for j = 1:size(angle_sum,2) 
        summary_data_delta(1,i,j) = prctile(all_medians(i,j,:),25);
        summary_data_delta(2,i,j) = prctile(all_medians(i,j,:),75);        
    end
end

summary_data=summary_data(2,:);
iqr_values = [];
plot_data_indiv = [];
plot_data_summary = [];
figure('Name', fig_name,'NumberTitle','off'); 
w = 1:2:6;
for s = 1:size(summary_data_median, 3)
    v = 0;   
    for t = 1:size(summary_data_median, 2)
        % median angle plot intact vs. inactivation
        subplot(size(summary_data, 2),2,w(s)); hold on 
        for u = 1:size(summary_data_median, 1)
            if u == 1
                scatter(u+v, summary_data_median(u, t, s), [], [0 0 0],"filled") 
                errorbar(u+v,  summary_data_median(u, t, s), summary_data_median(u, t, s) - summary_data_delta(1,u,t,s), summary_data_delta(2,u,t,s) - summary_data_median(u, t, s), 'LineWidth', 1, 'color', [0 0 0])
            elseif u == 2
                scatter(u+v, summary_data_median(u, t, s), [], [0.3010 0.7450 0.9330],"filled") 
                errorbar(u+v,  summary_data_median(u, t, s), summary_data_median(u, t, s) - summary_data_delta(1,u,t,s), summary_data_delta(2,u,t,s) - summary_data_median(u, t, s), 'LineWidth', 1, 'color', [0.3010 0.7450 0.9330])
            end
            plot_data_summary = [plot_data_summary;[u+v,  summary_data_median(u, t, s), summary_data_delta(1,u,t,s), summary_data_delta(2,u,t,s)]];
            temp_array=summary_data{1,t};
            scatter(u+v, temp_array(:,u), [], [.7 .7 .7],"filled")
            iqr_values = [iqr_values,summary_data_median(u, t, s) - summary_data_delta(1,u,t,s),summary_data_delta(2,u,t,s) - summary_data_median(u, t, s) ];
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

% signrank tests for the effects of laser, if needed
angle_sig_laser = zeros(3, 1);
angle_sig_position = zeros(2, 3);
for i = 1:size(angle_sig_laser, 1)
    for j = 1:size(angle_sig_laser, 2)
        [angle_sig_laser(i, j), ~, ~] = signrank(summary_data{i}(:,1), summary_data{i}(:,2), 'alpha', 0.05, 'tail', 'both');        
        for k = 1:size(angle_sig_position, 1)
            [angle_sig_position(k, i, j), ~, ~] = signrank(summary_data{i}(:,j),summary_data{2}(:,j), 'alpha', 0.05, 'tail', 'both');
        end
    end
end
end