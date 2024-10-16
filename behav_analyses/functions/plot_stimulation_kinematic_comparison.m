function [plot_data_summary,plot_data_indiv,shuffled_test_results,p_value_tailTwo_summary] = plot_stimulation_kinematics_comparison(data_set_cued,data_set_UNcued,ylabel_input);

load(data_set_cued)
data_cued = t_stats_combined;
load(data_set_uncued)
data_uncued = t_stats_combined;

gridSize=[5,3];
loop_num = 1000;
mm_pix = 0.06; % pixel to mm conversion
cued_raw=cell(gridSize);
ctrl_raw=cell(gridSize);
UNcued_raw=cell(gridSize);

% collect values for each of the three conditions
 for i=1:gridSize(1)
    for j=1:gridSize(2)
        temp_cued=cell(1,numel(data_cued));
        temp_ctrl=cell(1,numel(data_cued));
        temp_UNcued=cell(1,numel(data_cued));
        for k=1:numel(data_cued)
            for l = 1:2

                if l == 1
                    data_indiv = data_cued.data(k);
                elseif l == 2
                    data_indiv = data_uncued.data(k);
                end

                if ~isempty(data_indiv{i,j})
                    t_stats = data_indiv{i,j};
                    for p=1:numel(t_stats)
                        if strcmp(ylabel_input,'protLatency')
                            temp = t_stats(p).time_rel_cue;
                        elseif strcmp(ylabel_input,'duration')
                            temp = t_stats(p).dur;
                        elseif strcmp(ylabel_input,'InterLickInterval')
                            if p < numel(t_stats) && t_stats(p+1).laser_2D == 1 && t_stats(p+1).lick_index == lick_num+1 && isnan(t_stats(p).spout_contact)
                                temp = t_stats(p+1).time_rel_cue-t_stats(p).time_rel_cue;
                            end
                        elseif strcmp(ylabel_input,'pathlength')
                            temp = t_stats(p).pathlength_3D_t*0.06;
                        elseif strcmp(ylabel_input,'maxSpeed')                                
                            temp = max(t_stats(p).magspeed_t)*0.06*1000;
                        elseif strcmp(ylabel_input,'accPeakNum')
                            temp = numel(t_stats(p).accel_peaks_tot_t);
                        end                    

                        if l == 1 && t_stats(p).lick_index == lick_num && t_stats(p).laser_trial == 1 &&t_stats(p).laser_2D == 1 && t_stats(p).time_rel_cue>0
                            temp_cued{k}=[temp_cued{k},temp];
                        elseif l == 1 && t_stats(p).lick_index == lick_num && t_stats(p).laser_trial == 0 &&t_stats(p).laser_2D == 0 && t_stats(p).time_rel_cue>0 
                            temp_ctrl{k}=[temp_ctrl{k},temp];
                        elseif l == 2 && t_stats(p).lick_index == lick_num && t_stats(p).laser_trial == 1 &&t_stats(p).laser_2D == 1 && t_stats(p).time_rel_cue>0
                            temp_UNcued{k}=[temp_UNcued{k},temp];
                        end
                    end
                end
            end
        end
        cued_raw{i,j}=temp_cued;
        ctrl_raw{i,j}=temp_ctrl;
        UNcued_raw{i,j}=temp_UNcued;
    end
 end
 
% bin all sites together but not all mice
start_column=1; %first column to use
end_column=2; %last column to use
start_row=2; %first row to use
end_row=5; %last row to use
stats_cued=cell(1,numel(dirList_cued));
stats_ctrl=cell(1,numel(dirList_cued));
stats_UNcued=cell(1,numel(dirList_cued));
for i=start_row:end_row
    for j=start_column:end_column
        for k=1:numel(dirList_cued)
            stats_ctrl{k}=[stats_ctrl{k},ctrl_raw{i,j}{k}];
            stats_cued{k}=[stats_cued{k},cued_raw{i,j}{k}];
            stats_UNcued{k}=[stats_UNcued{k},UNcued_raw{i,j}{k}];
        end
    end
end

% bootstrap test for differences between the three conditions
data_all = vertcat(stats_ctrl,stats_cued,stats_UNcued);
median_BS = nan(loop_num, size(data_all,1));
shuffled_test_results = nan(loop_num,2,3);
for m=1:loop_num
    data_BS_xAnimal=datasample([1:size(data_all,2)],size(data_all,2));
    data_temp = data_all(:,data_BS_xAnimal);
    data_temp = cellfun(@(x) datasample(x,numel(x)),data_temp,'UniformOutput',false);
    data_temp = cellfun(@median,data_temp,'UniformOutput',false);
    median_BS(m,:) = median(cell2mat(data_temp),2)';
    shuffled_test_results(m,1,1) = median_BS(m,1) <= median_BS(m,2); % ctrl vs stim w/ cue, tailOne
    shuffled_test_results(m,2,1) = median_BS(m,1) >= median_BS(m,2); % ctrl vs stim w/ cue, tailTwo
    shuffled_test_results(m,1,2) = median_BS(m,1) <= median_BS(m,3); % ctrl vs stim w/o cue, tailOne
    shuffled_test_results(m,2,2) = median_BS(m,1) >= median_BS(m,3); % ctrl vs stim w/o cue, tailTwo
    shuffled_test_results(m,1,3) = median_BS(m,2) <= median_BS(m,3); % stim w/ cue vs stim w/o cue, tailOne
    shuffled_test_results(m,2,3) = median_BS(m,2) >= median_BS(m,3); % stim w/ cue vs stim w/o cue, tailTwo
end
shuffled_test_results = squeeze(mean(shuffled_test_results,1));

% get median and bootstrapped IQR
stats_cued=cellfun(@median,stats_cued);
stats_ctrl=cellfun(@median,stats_ctrl);
stats_UNcued=cellfun(@median,stats_UNcued);
median_actual = [median(stats_ctrl),median(stats_cued),median(stats_UNcued)];
liqr_BS = prctile(median_BS, 25,1);
hiqr_BS = prctile(median_BS, 75,1);

% make the plots
x_ticks=[0.25,0.75,1.25];
plot_data_indiv=[];
figure
errorbar(x_ticks,median_actual,median_actual-liqr_BS,hiqr_BS-median_actual, '.');
hold on
plot_data=horzcat(stats_ctrl',stats_cued',stats_UNcued');
for i=1:size(plot_data,2)
    scatter(ones(1,size(plot_data,1))*x_ticks(i),plot_data(:,i),'.');
    plot_data_indiv=[plot_data_indiv;plot_data(:,i)'];
end
xlim([0,1.5])
xticks(x_ticks)
xticklabels({'Control','Stim. w/ cue','Stim. w/o cue'})
xlabel('ML coordinates')
box off
plot_data_summary=[median_actual',liqr_BS',hiqr_BS'];

end