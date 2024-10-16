function [plot_summary,plot_indiv] = plot_stimulation_xsite_cumulProb(data_set,stim_type,animal_ID,start_column,end_column,start_row,end_row, fig_name)

% for this analysis, we only considered activations with the same 400ms durations
load(data_set)
gridSize=[5,3];
lick_num=1;

data_raw=cell(gridSize);
pooled_raw=[];
pooled_indiv=cell(1,numel(t_stats_combined));
for i=start_row:end_row
    for j=start_column:end_column
        temp_data=cell(1,numel(t_stats_combined));
        for k=1:numel(t_stats_combined)
            data_indiv = t_stats_combined(k).data;
            if sum(animal_ID==k) ~=1 || isempty(data_indiv{i,j})
                continue
            end
                       
            video_descriptor = t_stats_combined(k).video_descriptor{i,j};
            % for stim. without cues, number of trials is number of videos    
            if strcmp(stim_type,'uncued')
                t_stats = data_indiv{i,j};
                session_latency = nan(1,numel(t_stats_combined));
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
            pooled_indiv{k}=[pooled_indiv{k},session_latency];
        end                
        data_raw{i,j}=temp_data;
        pooled_raw=[pooled_raw,horzcat(temp_data{:})];
    end
end

% plot the cumulative probability curve
edges = 0:1:1000;
prot_vect = histc(pooled_raw,edges);
plot_vect(i,:) = cumsum(prot_vect)/numel(pooled_raw);
figure('Name', fig_name,'NumberTitle','off'); 
plot(plot_vect(i,:),'b');
xlim([1 400]);
ylim([-0.1 1]);
set(gcf,'position',[560 210 400 400])

pooled_indiv=pooled_indiv(animal_ID);
loop_num=1000;
% bootstrap to get IQR
protProb_BS = nan(loop_num,1);
for m=1:loop_num
    data_BS_xAnimal=datasample([1:numel(pooled_indiv)],numel(pooled_indiv));
    data_temp = pooled_indiv(data_BS_xAnimal); % get data from bootstrapped animals
    data_temp = cellfun(@(x) datasample(x,numel(x)),data_temp,'UniformOutput',false); % sample each animal's data
    data_temp = cellfun(@(x) numel(find(~isnan(x)))/numel(x), data_temp);
    protProb_BS(m) = median(data_temp, 'omitnan');
end
plot_indiv = cellfun(@(x) numel(find(~isnan(x)))/numel(x), pooled_indiv);
median_all = median(plot_indiv);
iqr_l = prctile(protProb_BS,25);
iqr_h = prctile(protProb_BS,75);
plot_summary = [median_all,iqr_l,iqr_h];

figure('Name', fig_name,'NumberTitle','off'); 
hold on
errorbar(1,median_all,median_all-iqr_l,iqr_h-median_all,'.');
scatter(ones(1,numel(plot_indiv)),plot_indiv,120,'.');
ylim([0,1.1])
ylabel('Protrusion Probability')
box off
end