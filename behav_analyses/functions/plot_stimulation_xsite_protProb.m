function [plot_data_summary,plot_data_indiv,p_kw,c_kw] = plot_stimulation_xsite_protProb(data_set,stim_type,start_column,end_column,start_row,end_row, fig_name)

load(data_set)
gridSize=[5,3];
lick_num=1;
loop_num=1000;

data_raw=cell(gridSize);
for i=1:gridSize(1)
    for j=1:gridSize(2)
        temp_data=cell(1,numel(t_stats_combined));
        for k=1:numel(t_stats_combined)
            data_indiv = t_stats_combined(k).data;
            if isempty(data_indiv{i,j})
                continue
            end
                       
            video_descriptor = t_stats_combined(k).video_descriptor{i,j};
            % for stim. without cues, number of trials is number of videos   
            if strcmp(stim_type,'uncued')
                t_stats = data_indiv{i,j};
                session_latency = nan(1,numel(video_descriptor));
                for p=1:max([t_stats.trial_num])
                    l_index = find(([t_stats.trial_num]==p)&([t_stats.time_rel_cue]>0)&([t_stats.laser_2D]==1));
                    if isempty(l_index)
                        continue
                    end
                    l_stats = t_stats(l_index);
                    if isnan(find([l_stats.lick_index] == lick_num))
                        continue
                    end
                    session_latency(p)=l_stats([l_stats.lick_index] == lick_num).time_rel_cue;
                end
            % for stim. with cues or control trials, number of trials is number of videos tagged with and without L respectively
            elseif ~strcmp(stim_type,'uncued')               
                laser_trial_num = 0;
                for videoNum = 1:numel(video_descriptor)
                    descriptor = video_descriptor{videoNum};
                    if (strcmp(stim_type,'cued') && descriptor(end) == 'L') || (strcmp(stim_type,'ctrl') && descriptor(end) ~= 'L')
                        laser_trial_num = laser_trial_num + 1;
                    end
                end
                
                t_stats = data_indiv{i,j};               
                session_latency = nan(1,laser_trial_num);
                counter = 1;
                    for p=1:max([t_stats.trial_num])                        
                        if strcmp(stim_type,'cued')
                            l_index = find(([t_stats.trial_num]==p)&([t_stats.time_rel_cue]>0)&([t_stats.laser_2D]==1));
                        elseif strcmp(stim_type,'ctrl')
                            l_index = find(([t_stats.trial_num]==p)&([t_stats.time_rel_cue]>0)&([t_stats.laser_trial]==0));
                        end
                        if isempty(l_index)
                            continue
                        end
                        l_stats = t_stats(l_index);
                        if isempty(find([l_stats.lick_index] == lick_num))
                            continue
                        end
                        session_latency(counter)=l_stats([l_stats.lick_index] == lick_num).time_rel_cue;
                        counter = counter + 1;
                    end
            end            
            temp_data{k}=session_latency;
        end               
        data_raw{i,j}=temp_data;
        all_medians=nan(1,loop_num);
        temp_data=temp_data(~cellfun('isempty', temp_data));
    end
end

% combine data with the same AP for each mouse
AP_data=cell(1,gridSize(1));
AP_data_raw=cell(1,gridSize(1));
for i=1:gridSize(1)
    temp_AP_data=cell(1,numel(t_stats_combined));
    for j=1:numel(t_stats_combined)
        temp_data=[];
        for k=start_column:end_column
            temp_data=[temp_data,data_raw{i,k}{j}];
        end
        temp_AP_data{j}=temp_data;
    end
    AP_data_raw{i}=temp_AP_data;
    
    for q=1:numel(temp_AP_data)
        AP_data{i}=[AP_data{i},numel(find(~isnan(temp_AP_data{q})))/numel(temp_AP_data{q})];
    end
end

% if you want to pool data across sites, the following code gives you the pooled median and bootstrapped IQR of invidual mice as well as acoss mice
pooled_data_raw = vertcat(AP_data_raw{start_row:end_row});
for i = 1:size(pooled_data_raw,2)
    pooled_data_raw{1, i} = horzcat(pooled_data_raw{:,i});   
end
pooled_data_raw = pooled_data_raw(1, :);
pooled_median_BS = nan(loop_num,1);
for m=1:loop_num
    data_BS_xAnimal=datasample([1:numel(pooled_data_raw)],numel(pooled_data_raw));
    data_temp = pooled_data_raw(data_BS_xAnimal);
    data_temp = cellfun(@(x) datasample(x,numel(x)),data_temp,'UniformOutput',false);
    data_temp = cellfun(@(x) numel(find(~isnan(x)))/numel(x), data_temp);
    pooled_median_BS(m) = median(data_temp, 'omitnan');
end
pooled_individual = cellfun(@(x) numel(find(~isnan(x)))/numel(x), pooled_data_raw);

% bootstrap to get IQR for each site
AP_median_BS = nan(loop_num,numel(AP_data_raw));
for m=1:loop_num
    data_BS_xAnimal=datasample([1:numel(AP_data_raw{1,1})],numel(AP_data_raw{1,1}));
    for n=1:numel(AP_data_raw)
        data_temp = AP_data_raw{n}(data_BS_xAnimal); % get data from bs animals for a stim site
        data_temp = cellfun(@(x) datasample(x,numel(x)),data_temp,'UniformOutput',false); % sample each animal's data
        data_temp = cellfun(@(x) numel(find(~isnan(x)))/numel(x), data_temp);
        AP_median_BS(m,n) = median(data_temp, 'omitnan');
    end
end
AP_median = cellfun(@(x) median(x,'omitnan'), AP_data,'UniformOutput',false);
AP_iqr_h=prctile(AP_median_BS,75,1);
AP_iqr_l=prctile(AP_median_BS,25,1);
AP_data=AP_data(start_row:end_row);
AP_median=AP_median(start_row:end_row);
AP_iqr_l=AP_iqr_l(start_row:end_row);
AP_iqr_h=AP_iqr_h(start_row:end_row);

% plot data across sites
figure('Name', fig_name,'NumberTitle','off'); 
xval= [3.1 3.4 3.7 4];
hold on
for i=1:size(AP_data,2)
    errorbar(xval(i),AP_median{i},AP_median{i}-AP_iqr_l(i),AP_iqr_h(i)-AP_median{i},'.');
    scatter(ones(1,numel(AP_data{1,i}))*xval(i),AP_data{1,i},120,'.');
end
xlim([min(xval)-0.15,max(xval)+0.15])
ylim([0,1.1])
xticks(xval)
xticklabels({'-3.1AP','-3.4AP','-3.7AP','-3.9AP'})
yticks([0,0.5,1])
xlabel('AP coordinates')
ylabel('Protrusion Probability')
box off

for i = 1:numel(AP_median)
    plot_data_summary(i,1) = AP_median{i};
    plot_data_summary(i,2) = AP_iqr_l(i);
    plot_data_summary(i,3) = AP_iqr_h(i);
    plot_data_indiv(i,:) = AP_data{i};
end

% correlation test
test_data=AP_data;
stats_data=[];
site_data=[];
for i=1:size(test_data,2)
    temp_data=test_data{i};
    nanIND=isnan(temp_data);
    stats_data=[stats_data;temp_data(~nanIND)'];
    site_data=[site_data;i*ones(numel(temp_data(~nanIND)),1)];
end
[r,p] = corr(site_data,stats_data, 'Type', 'Spearman');

%kruskal wallis test;
stats_data=[];
for i=1:size(test_data,2)
    stats_data=[stats_data,test_data{1,i}'];
end
[p_kw,tbl,stats] = kruskalwallis(stats_data);
c_kw = multcompare(stats,'Display','off');

end