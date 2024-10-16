function plot_angle_correlation_example(t_stats_paths,lick_1,lick_2,indiv_trial_type, fig_name)

if ~iscell(t_stats_paths)
    % If only a single path is passed in, wrap it in a cell array
    t_stats_paths = {t_stats_paths};
end

% Find mouse IDs for all sessions in t_stats_paths
mouse_ID_session = cell(numel(t_stats_paths), 1);
for i = 1:numel(t_stats_paths)
    mouse_ID_session{i} = t_stats_paths{i}(regexp(t_stats_paths{i}, '(\w)*(?=\\Masks)'):regexp(t_stats_paths{i}, '\w(?=\\Masks)'));
end
clear i

% find unique mice
unique_mouse_ID = unique(mouse_ID_session);
lick_1_fakeout = {};
mouse_ID = {};
ml_error_ang = {};
prot_angle_change = {};

for i = 1:numel(unique_mouse_ID)
    
    % find indices in mouse_ID_session corresponding to a certain mouse
    mouse_session_ind = find(strcmp(mouse_ID_session, unique_mouse_ID{i}));
    
    % concatenate sessions from same mice
    t_stats_all = [];
    fiducial_xy_all = [];
    for j = 1:sum(strcmp(mouse_ID_session, unique_mouse_ID{i}))
        load(t_stats_paths{mouse_session_ind(j)})
        t_stats_all = [t_stats_all, t_stats];
        fiducial_xy = [(-1)*unique(t_stats(1).tip_xf-t_stats(1).tip_x)+(400 - 240), (-1)*unique(t_stats(i).tip_yf-t_stats(i).tip_y)];
        fiducial_xy_all = [fiducial_xy_all, repmat(fiducial_xy(mouse_session_ind(j), :), numel(t_stats), 1)];
    end
    clear j
          
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

ml_error_temp = vertcat(ml_error_ang{:});
prot_angle_change_temp = vertcat(prot_angle_change{:});
lick_1_fakeout_temp = vertcat(lick_1_fakeout{:});

% fit robust regression line
[b, stats] = robustfit(ml_error_temp, prot_angle_change_temp);

% plot for all conditions combined
marker_size = 50;

figure('Name', fig_name,'NumberTitle','off');
scatter(ml_error_temp(lick_1_fakeout_temp == 1), prot_angle_change_temp(lick_1_fakeout_temp == 1), marker_size, [247/255 148/255 30/255]);
hold on;
scatter(ml_error_temp(lick_1_fakeout_temp == 2), prot_angle_change_temp(lick_1_fakeout_temp == 2), marker_size, [0 0 0]);
scatter(ml_error_temp(lick_1_fakeout_temp == 3), prot_angle_change_temp(lick_1_fakeout_temp == 3), marker_size, [127/255 63/255 152/255]);
y1 = b(1) + b(2)*ml_error_temp;
plot(ml_error_temp, y1);
xlabel('ML Error')
ylabel('Change in Protrusion Angle from L3 to L2')

[r, p] = corrcoef(ml_error_temp, prot_angle_change_temp);

% plot for each condition
colors = {[247/255 148/255 30/255],[0 0 0], [127/255 63/255 152/255]};
marker_size = 50;
for i = 1:numel(indiv_trial_type)
    lick_num = indiv_trial_type(i);
    [b, stats] = robustfit(ml_error_temp(lick_1_fakeout_temp == lick_num), prot_angle_change_temp(lick_1_fakeout_temp == lick_num));

    figure('Name', fig_name,'NumberTitle','off'); 
    scatter(ml_error_temp(lick_1_fakeout_temp == lick_num), prot_angle_change_temp(lick_1_fakeout_temp == lick_num), marker_size, colors{lick_num});
    hold on;
    y1 = b(1) + b(2)*ml_error_temp;
    plot(ml_error_temp, y1);
    [r, p] = corrcoef(ml_error_temp(lick_1_fakeout_temp == lick_num), prot_angle_change_temp(lick_1_fakeout_temp == lick_num));
end
end