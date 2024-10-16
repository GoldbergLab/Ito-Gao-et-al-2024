function [shuffle_test_lesion_combined_tailOne,shuffle_test_lesion_combined_tailTwo,plot_data_summary_combined,plot_data_indiv_combined] = plot_xsession_comparison_initiation(data_set_one,data_set_two, fig_name)

load(data_set_one)
animal_t_stats_paths_preLesion = t_stats_combined;
load(data_set_two)
animal_t_stats_paths_postLesion = t_stats_combined;

initiation_all = cell(numel(animal_t_stats_paths_preLesion),2);
CSM_all = cell(numel(animal_t_stats_paths_preLesion),2); % this script collects CSM data but currently does not analyze it
for p = 1:numel(animal_t_stats_paths_preLesion)
    for q = 1:2
        if q == 1
            t_stats = animal_t_stats_paths_preLesion(p).data;
        elseif q == 2
            t_stats = animal_t_stats_paths_postLesion(p).data;
        end
        temp = zeros(1,max([t_stats.trial_num]));
        temp_CSM = nan(1,max([t_stats.trial_num]));
        for i = 1:numel(t_stats)
            if t_stats(i).lick_index == 1 && t_stats(i).time_rel_cue>0 
                temp(t_stats(i).trial_num) = 1;
                if t_stats(i).CSM_dur > 5
                    temp_CSM(t_stats(i).trial_num) = 1;
                else
                    temp_CSM(t_stats(i).trial_num) = 0;
                end
            end
        end
        lastTrial = find(temp,1,'last');
        temp = temp(1:lastTrial);
        initiation_all{p,q} = temp;
        CSM_all{p,q} = temp_CSM(~isnan(temp_CSM));
    end
end

initiation_by_mouse = cellfun(@mean, initiation_all);

numShuffles = 1000;
shuf_initiation = zeros([numShuffles size(initiation_all)]);
for i = 1:numShuffles

    % randomly select mice to create a new pseudopopulation
    shuf_animal_ind = randi(size(initiation_all, 1), 1, size(initiation_all, 1));
    shuf_initiation_animal = initiation_all(shuf_animal_ind, :);
    
    % randomly select trials from each mouse/condition to get a new pseudosession
    for j = 1:size(shuf_initiation_animal, 1)
        for k = 1:size(shuf_initiation_animal, 2)
            trial_num_temp = cellfun(@numel, squeeze(shuf_initiation_animal(j, k, :, :)), 'UniformOutput', false);
            
            if sum([trial_num_temp{:}]) == 0
                shuf_initiation(i, j, k) = [NaN;NaN;NaN];
            else
                trial_num_ind_temp = cellfun(@(x) randi(x, 1, x), trial_num_temp, 'UniformOutput', false);
                shuf_initiation_temp = cellfun(@(y, z) y(z), squeeze(shuf_initiation_animal(j, k, :, :)), trial_num_ind_temp, 'UniformOutput', false);
                shuf_initiation(i, j, k) = cellfun(@mean, shuf_initiation_temp);
            end
        end
    end
end

shuf_initiation_all = squeeze(median(shuf_initiation, 2));

shuffle_test_lesion_combined_tailOne = nan(numShuffles,1);
shuffle_test_lesion_combined_tailTwo = nan(numShuffles,1);
for i = 1:numShuffles
    shuffle_test_lesion_combined_tailOne(i) = shuf_initiation_all(i,1) - shuf_initiation_all(i,2) <= 0;
    shuffle_test_lesion_combined_tailTwo(i) = shuf_initiation_all(i,1) - shuf_initiation_all(i,2) >= 0;
end
shuffle_test_lesion_combined_tailOne = mean(shuffle_test_lesion_combined_tailOne,1);
shuffle_test_lesion_combined_tailTwo = mean(shuffle_test_lesion_combined_tailTwo,1);

figure('Name', fig_name,'NumberTitle','off'); 
hold on
for u = 1:size(initiation_by_mouse,2)
    if u == 1
        scatter(1.2*ones(1,size(initiation_by_mouse,1)), initiation_by_mouse(:,u), [], [0 0 0],"filled") 
        errorbar(1.2, median(initiation_by_mouse(:,u)), median(initiation_by_mouse(:,u))-prctile(shuf_initiation_all(:,u),25), prctile(shuf_initiation_all(:,u),75)-median(initiation_by_mouse(:,u)), 'LineWidth', 1, 'color', [0 0 0])
    elseif u == 2
        scatter(1.8*ones(1,size(initiation_by_mouse,1)), initiation_by_mouse(:,u), [], [0.3010 0.7450 0.9330],"filled") 
        errorbar(1.8, median(initiation_by_mouse(:,u)), median(initiation_by_mouse(:,u))-prctile(shuf_initiation_all(:,u),25), prctile(shuf_initiation_all(:,u),75)-median(initiation_by_mouse(:,u)), 'LineWidth', 1, 'color', [0.3010 0.7450 0.9330])
    end
end
xlim([1 2])
ylim([-0.05 1.05])
xticks([1.2,1.8])
xticklabels({'Pre-lesion', 'Post-lesion'})
ylabel('Lick Initiation Rate')

plot_data_summary_combined = [median(initiation_by_mouse(:,1)), prctile(shuf_initiation_all(:,1),25), prctile(shuf_initiation_all(:,1),75);...
                              median(initiation_by_mouse(:,2)), prctile(shuf_initiation_all(:,2),25), prctile(shuf_initiation_all(:,2),75)];
plot_data_indiv_combined = [initiation_by_mouse(:,1);initiation_by_mouse(:,2)];
end