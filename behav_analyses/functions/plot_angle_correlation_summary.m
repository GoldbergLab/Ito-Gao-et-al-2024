function [plot_data_summary,plot_data_indiv,p_value_tailOne_summary,p_value_tailTwo_summary] = plot_angle_correlation_summary(data_set, lick_1, lick_2, indiv_trial_type, fig_name)

load(data_set)
% Below is to loop across all sessions and summary data across mice
% mouse_ID_session = cell(numel(t_stats_combined), 1);
% for i = 1:numel(t_stats_combined)
%     mouse_ID_session{i} = t_stats_paths{i}(regexp(t_stats_paths{i}, '(\w)*(?=\\Masks)'):regexp(t_stats_paths{i}, '\w(?=\\Masks)'));
% end
% clear i

% find unique mice
%unique_mouse_ID = unique(mouse_ID_session);
lick_1_fakeout = {};
mouse_ID = {};
ml_error = {};
ml_error_ang = {};
prot_angle_change = {};

for i = 1:numel(t_stats_combined)
    
    % find indices in mouse_ID_session corresponding to a certain mouse
    %mouse_session_ind = find(strcmp(mouse_ID_session, unique_mouse_ID{i}));
    
    % concatenate sessions from same mice
    % t_stats_all = [];
    % fiducial_xy_all = [];
    % for j = 1:sum(strcmp(mouse_ID_session, unique_mouse_ID{i}))
    %     load(t_stats_paths{mouse_session_ind(j)})
    %     t_stats_all = [t_stats_all, t_stats];
    %     fiducial_xy = [(-1)*unique(t_stats(1).tip_xf-t_stats(1).tip_x)+(400 - 240), (-1)*unique(t_stats(i).tip_yf-t_stats(i).tip_y)];
    %     fiducial_xy_all = [fiducial_xy_all, repmat(fiducial_xy, numel(t_stats), 1)];
    % end
    % clear j

    %load(t_stats_combined{i}.data)
    t_stats = t_stats_combined(i).data;
    t_stats_all = t_stats;
    fiducial_xy = [(-1)*unique(t_stats(1).tip_xf-t_stats(1).tip_x)+(400 - 240), (-1)*unique(t_stats(i).tip_yf-t_stats(i).tip_y)];
    fiducial_xy_all = repmat(fiducial_xy, numel(t_stats), 1);
          
    % Find trials in which we want to analyze
    contact_centroid_exist = unique([t_stats_all([t_stats_all.lick_index_contact1] == lick_1 & cellfun(@(x) numel(x) ~= 1, {t_stats_all.contact_centroid}) & ~isnan([t_stats_all.spout_contact])).trial_num]);
    lick_2_exist = unique([t_stats_all([t_stats_all.lick_index_contact1] == lick_2).trial_num]);
    lick_1_no_double_tap = unique([t_stats_all([t_stats_all.lick_index_contact1] == 1 & isnan([t_stats_all.spout_contact2])).trial_num]);
    no_laser = unique([t_stats_all([t_stats_all.laser_trial] == 0).trial_num]);    
    cont_first_ind = [t_stats_all(([t_stats_all.spout_contact] - [t_stats_all.time_rel_cue] > 0) & [t_stats_all.lick_index_contact1] == lick_1).trial_num];
    lick_1_contact_centroid_exist = unique([t_stats_all([t_stats_all.lick_index_contact1] == lick_1 & cellfun(@(x) ~isnan(x(1)), {t_stats_all.contact_centroid}) & ~isnan([t_stats_all.spout_contact])).trial_num]);
    trial_nums = intersect(lick_1_contact_centroid_exist, intersect(cont_first_ind, intersect(no_laser, intersect(lick_1_no_double_tap, intersect(contact_centroid_exist, lick_2_exist)))));

    % Filter t_stats accordingly
    t_stats_temp = t_stats_all(ismember([t_stats_all.trial_num], trial_nums));
    fiducial_xy_temp = fiducial_xy_all(ismember([t_stats_all.trial_num], trial_nums), :);
    fiducial_xy_temp(:, 1) = fiducial_xy_temp(:, 1) - (400 - 240);

    % Get fakeout trial type
    lick_1_fakeout_temp = [t_stats_temp([t_stats_temp.lick_index_contact1] == lick_1).fakeout_trial]';

    % Get L2 contact centroid
    lick_1_cont_centroid_temp = cellfun(@(x) x(1, 1:2), {t_stats_temp([t_stats_temp.lick_index_contact1] == lick_1).contact_centroid}, 'UniformOutput', false);
    lick_1_cont_centroid_temp = cell2mat(lick_1_cont_centroid_temp(:));
    lick_1_cont_centroid_temp(:, 1) = lick_1_cont_centroid_temp(:, 1) - fiducial_xy_temp([t_stats_temp.lick_index_contact1] == lick_1, 1);
    lick_1_cont_centroid_temp(:, 2) = lick_1_cont_centroid_temp(:, 2) - fiducial_xy_temp([t_stats_temp.lick_index_contact1] == lick_1, 2);
    
    % Get L2 lateral displacement
    lick_1_tip_x = {t_stats_temp([t_stats_temp.lick_index_contact1] == lick_1).tip_x}';
    lick_1_tip_y = {t_stats_temp([t_stats_temp.lick_index_contact1] == lick_1).tip_y}';
    lick_1_time_rel_cue = {t_stats_temp([t_stats_temp.lick_index_contact1] == lick_1).time_rel_cue}';
    lick_1_spout_contact = {t_stats_temp([t_stats_temp.lick_index_contact1] == lick_1).spout_contact}';
    lick_1_cont_ind = num2cell(cellfun(@(x, y) x - y, lick_1_spout_contact, lick_1_time_rel_cue));
    lick_1_tip_xy_temp = [cellfun(@(x, y) x(y), lick_1_tip_x, lick_1_cont_ind) cellfun(@(x, y) x(y), lick_1_tip_y, lick_1_cont_ind)];
    lick_1_tip_xy_temp(:, 1) = lick_1_tip_xy_temp(:, 1) - fiducial_xy_temp(1, 1);
    lick_1_tip_xy_temp(:, 2) = lick_1_tip_xy_temp(:, 2) - fiducial_xy_temp(1, 2);
    
    % Get L3 lateral displacement
    lick_2_tip_x = {t_stats_temp([t_stats_temp.lick_index_contact1] == lick_2).tip_x}';
    lick_2_tip_y = {t_stats_temp([t_stats_temp.lick_index_contact1] == lick_2).tip_y}';
    lick_2_prot_ind = {t_stats_temp([t_stats_temp.lick_index_contact1] == lick_2).prot_ind}';
    lick_2_tip_xy_temp = [cellfun(@(x, y) x(y), lick_2_tip_x, lick_2_prot_ind) cellfun(@(x, y) x(y), lick_2_tip_y, lick_2_prot_ind)];
    lick_2_tip_xy_temp(:, 1) = lick_2_tip_xy_temp(:, 1) - fiducial_xy_temp(1, 1);
    lick_2_tip_xy_temp(:, 2) = lick_2_tip_xy_temp(:, 2) - fiducial_xy_temp(1, 2);
    
    % Create mouse_ID
    mouse_ID_temp = repmat(i, numel(trial_nums), 1);
    
    lick_1_fakeout{i} = lick_1_fakeout_temp;
    mouse_ID{i} = mouse_ID_temp;
        
    % error angle as an angle between two vectors
    x1=lick_1_cont_centroid_temp(:, 1);
    y1=lick_1_cont_centroid_temp(:, 2); 
    x2=lick_1_tip_xy_temp(:, 1);
    y2=lick_1_tip_xy_temp(:, 2);
    ml_error_ang{i} = atan2d(x1.*y2-y1.*x2,x1.*x2+y1.*y2);
    
    % calculate protrusion angle change from lick_1 to lick_2
    prot_angle_change{i} = asind(lick_2_tip_xy_temp(:, 1)./sqrt(lick_2_tip_xy_temp(:, 1).^2 + lick_2_tip_xy_temp(:, 2).^2)) - asind(lick_1_tip_xy_temp(:, 1)./sqrt(lick_1_tip_xy_temp(:, 1).^2 + lick_1_tip_xy_temp(:, 2).^2));
end

% The code below shuffles and plots for across conditions and mice.
r_value = [];
p_value = [];
for i = 1:numel(ml_error_ang)
    ml_error_temp = vertcat(ml_error_ang{i});
    prot_angle_change_temp = vertcat(prot_angle_change{i});
    
    [r, p] = corrcoef(ml_error_temp, prot_angle_change_temp);
    
    r_value(i) = r(1, 2);
    p_value(i) = p(1, 2);    
end

% shuffle to create a null distribution 
numShuffles = 1000;
r_value_null = nan(numShuffles, 1);
for j = 1:numShuffles    
    prot_angle_change_shuf = cellfun(@(x) x(randperm(numel(x))), prot_angle_change, 'UniformOutput', false); %JG changed this line
    
    r_value_shuf_temp = [];
    for k = 1:numel(ml_error_ang)
        [r_temp, ~] = corrcoef(vertcat(ml_error_ang{k}), vertcat(prot_angle_change_shuf{k}));
        r_value_shuf_temp(k) = r_temp(1, 2);
    end
    r_value_null(j) = median(r_value_shuf_temp);
end

% shuffle to get IQR for mean correlation coefficient
numShuffles = 1000;
r_value_shuf = nan(numShuffles, 1);
for j = 1:numShuffles
    ml_error_ang_shuf_id = randi(numel(ml_error_ang), numel(ml_error_ang), 1);
    ml_error_ang_shuf = ml_error_ang(ml_error_ang_shuf_id);
    prot_angle_change_shuf = prot_angle_change(ml_error_ang_shuf_id);
    
    r_value_shuf_temp = [];
    for k = 1:numel(ml_error_ang_shuf)
        [r_temp, ~] = corrcoef(vertcat(ml_error_ang_shuf{k}), vertcat(prot_angle_change_shuf{k}));
        r_value_shuf_temp(k) = r_temp(1, 2);
    end
    
    r_value_shuf(j) = median(r_value_shuf_temp);    
end

r_value_iqr = prctile(r_value_shuf, [25 75]);
r_value_median = median(r_value);

r_null_iqr = prctile(r_value_null, [25 75]);
r_null_median = median(r_value_null);

plot_data_indiv = [];
plot_data_summary = [];
p_value_tailOne_summary = [];
p_value_tailTwo_summary = [];
figure('Name', fig_name,'NumberTitle','off'); 
scatter(1, r_value_median);
hold on;
errorbar(1, r_value_median, r_value_median - r_value_iqr(1), r_value_iqr(2) - r_value_median,'.');
scatter(3, r_null_median,'.');
scatter(ones(1,numel(r_value)),r_value,'.')
errorbar(3, r_null_median, r_null_median - r_null_iqr(1), r_null_iqr(2) - r_null_median,'.');
plot_data_summary = [r_value_median,r_value_iqr(1),r_value_iqr(2);...
                     r_null_median,r_null_iqr(1),r_null_iqr(2)];
plot_data_indiv = [r_value];
p_value_tailOne_summary = [numel(find(r_value_null<=r_value_median))/numel(r_value_null)];
p_value_tailTwo_summary = [numel(find(r_value_null>=r_value_median))/numel(r_value_null)];
title('Correlation coefficient across all mice');
xlim([0 4]);
ylim([-0.1 1]);
yticks([0 1]);
xticks([1 3]);

% The code below shuffles and plots for within condition acrosss mice.
% indiv_trial_type = [1,2,3]; % which conditions to plot individually  bmk27 made this a function parameter 10/16/2024

r_value = [];
p_value = [];
for i = 1:numel(ml_error_ang)
    ml_error_temp = vertcat(ml_error_ang{i});
    prot_angle_change_temp = vertcat(prot_angle_change{i});
    lick_1_fakeout_temp = vertcat(lick_1_fakeout{i});
    
    for j = 1:3
        [r, p] = corrcoef(ml_error_temp(lick_1_fakeout_temp == j), prot_angle_change_temp(lick_1_fakeout_temp == j));
        r_value(i, j) = r(1, 2);
        p_value(i, j) = p(1, 2);    
    end
end

% shuffle to create a null distribution 
numShuffles = 1000;
r_value_null = nan(numShuffles, 3);
r_value_shuf_temp = [];
for j = 1:numShuffles    
    for k = 1:numel(ml_error_ang)
        for m = 1:3
            ml_error_temp = ml_error_ang{k}(lick_1_fakeout{k} == m);
            prot_angle_change_temp = prot_angle_change{k}(lick_1_fakeout{k} == m);
            
            prot_angle_change_shuf = randperm(numel(prot_angle_change_temp));
            [r_temp, ~] = corrcoef(ml_error_temp, prot_angle_change_shuf);
            r_value_shuf_temp(k, m) = r_temp(1, 2);
        end
    end
    
    r_value_null(j, :) = median(r_value_shuf_temp, 1);
end

% shuffle to get IQR for mean correlation coefficient
numShuffles = 1000;
r_value_shuf_temp = [];
r_value_shuf = nan(numShuffles, 3);
for j = 1:numShuffles
    ml_error_ang_shuf_id = randi(numel(ml_error_ang), numel(ml_error_ang), 1);
    ml_error_ang_shuf = ml_error_ang(ml_error_ang_shuf_id);
    prot_angle_change_shuf = prot_angle_change(ml_error_ang_shuf_id);
    lick_1_fakeout_shuf = lick_1_fakeout(ml_error_ang_shuf_id);
    
    for k = 1:numel(ml_error_ang_shuf)
        for m = 1:3
            [r_temp, ~] = corrcoef(vertcat(ml_error_ang_shuf{k}(lick_1_fakeout_shuf{k} == m)), vertcat(prot_angle_change_shuf{k}(lick_1_fakeout_shuf{k} == m)));
            r_value_shuf_temp(k, m) = r_temp(1, 2);
        end
    end
    
    r_value_shuf(j, :) = median(r_value_shuf_temp, 1);
    
end

r_value_iqr = prctile(r_value_shuf, [25 75], 1);
r_value_median = median(r_value);
r_null_iqr = prctile(r_value_null, [25 75], 1);
r_null_median = median(r_value_null);

titles = {'Corr-coeff across left trials','Corr-coeff across center trials', 'Corr-coeff across right trials'};
for i = 1:numel(indiv_trial_type)
    fakeout = indiv_trial_type(i);
    figure('Name', fig_name,'NumberTitle','off'); 
    scatter(1, r_value_median(fakeout));
    hold on;
    errorbar(1, r_value_median(fakeout), r_value_median(fakeout) - r_value_iqr(1, fakeout), r_value_iqr(2, fakeout) - r_value_median(fakeout),'.');
    scatter(ones(1,numel(r_value(:,fakeout))),r_value(:,fakeout),'.')
    scatter(3, r_null_median(fakeout),'.');
    errorbar(3, r_null_median(fakeout), r_null_median(fakeout) - r_null_iqr(1, fakeout), r_null_iqr(2, fakeout) - r_null_median(fakeout),'.');
    title(titles{i})
    xlim([0 4]);
    ylim([-0.1 1]);
    yticks([0 1]);
    xticks([1 3]);
    plot_data_summary = [plot_data_summary;r_value_median(fakeout),r_value_iqr(1, fakeout),...
                        r_value_iqr(2, fakeout);r_null_median(fakeout),r_null_iqr(1, fakeout),...
                        r_null_iqr(2, fakeout)];
    plot_data_indiv = [plot_data_indiv;r_value(:,fakeout)'];
    p_value_tailOne_summary = [p_value_tailOne_summary,numel(find(r_value_null(fakeout,:)>=r_value_median(fakeout)))/numel(r_value_null(fakeout,:))];
    p_value_tailTwo_summary = [p_value_tailTwo_summary,numel(find(r_value_null(fakeout,:)<=r_value_median(fakeout)))/numel(r_value_null(fakeout,:))];

end
end