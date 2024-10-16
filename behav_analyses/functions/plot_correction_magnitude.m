function [plot_data_summary,plot_data_indiv,shuffle_test_laser_tailOne,shuffle_test_laser_tailTwo,shuffle_test_brain_tailOne,shuffle_test_brain_tailTwo] = plot_correction_magnitude(angle_sum,ylim_val, fig_name)
summary_data = angle_sum;
reverse_index=[size(summary_data,2):-1:1];

% organize data as ipsi/center/right, instead of left/center/reight, reverse +/- sign if the inactivation is on the right hemisphere
for i = (size(summary_data,1)/2) + 1 : size(summary_data,1)
    data_temp = {summary_data{i,:}};
    for j = 1:size(summary_data,2)
        summary_data{i,reverse_index(j)} = cellfun(@(x) x*(-1), data_temp{1,j},'un',0);
    end
end

% Combine data from the same mouse, note that some mice may only have unilateral inactivation on one hemisphere. The following code takes care of that
summary_data_median = [];
for i = 1:(size(summary_data,1)/2)
    for j = 1:size(summary_data,2)
        ipsi_data = summary_data{i,j};
        contra_data = summary_data{i+(size(summary_data,1)/2),j};
        for k = 1:size(summary_data{i,j},2)
            for l = 1:min(size(ipsi_data,1), size(contra_data,1))
                summary_data{i,j}{l,k} = vertcat(ipsi_data{l,k}, contra_data{l,k});
            end
        end
        if size(ipsi_data,1) < size(contra_data,1)
            summary_data{i,j} = vertcat(summary_data{i,j},contra_data(l+1:end,:));
        end
    end
end

% shuffle to find the IQR and perform stats tests
numShuffles = 1000;
all_median = nan(size(summary_data,1)/2,2,2,numShuffles); % brain x side x laser
shuffle_test_laser_tailOne = nan(size(summary_data,1)/2,2,numShuffles); % brain region x side
shuffle_test_laser_tailTwo = nan(size(summary_data,1)/2,2,numShuffles);
shuffle_test_brain_tailOne = nan(nchoosek(size(summary_data,1)/2,2),2,numShuffles); % combination of brain regions x side
shuffle_test_brain_tailTwo = nan(nchoosek(size(summary_data,1)/2,2),2,numShuffles);

for i = 1:numShuffles
    for j = 1:(size(summary_data,1)/2)
        %shuf_animal_ind = randi(size(summary_data{j}, 1), 1, size(summary_data{j}, 1)); % resample animals with replacement
        shuf_animal_ind = datasample([1:size(summary_data{j}, 1)],size(summary_data{j}, 1));
        
        for k = 1:size(summary_data{j}, 2)
            shuf_angle_animal = horzcat(summary_data{j,1}(shuf_animal_ind, k),summary_data{j,2}(shuf_animal_ind, k),summary_data{j,3}(shuf_animal_ind, k));
            for m = 1:size(shuf_angle_animal,1)
                for n = 1:size(shuf_angle_animal,2)
                    data_temp = shuf_angle_animal{m,n};
                    shuf_angle_animal{m,n} = datasample(data_temp,numel(data_temp)); % resample trials within each condition with replacement
                end
            end
            all_median(j,1,k,i) = median(cellfun(@median, shuf_angle_animal(:,2)) - cellfun(@median, shuf_angle_animal(:,1)));
            all_median(j,2,k,i) = median(cellfun(@median, shuf_angle_animal(:,3)) - cellfun(@median, shuf_angle_animal(:,2)));
        end
    end

    for j = 1:size(shuffle_test_laser_tailOne,2)
        for k = 1:size(shuffle_test_laser_tailOne,1)
            % laser off vs laser on
            shuffle_test_laser_tailOne(k,j,i) = all_median(k,j,1,i) - all_median(k,j,2,i) <= 0;
            shuffle_test_laser_tailTwo(k,j,i) = all_median(k,j,1,i) - all_median(k,j,2,i) >= 0;
        end
        % brain regions vs each other, for laser on only
        shuffle_test_brain_tailOne(1,j,i) = all_median(1,j,2,i) - all_median(2,j,2,i) <= 0;
        shuffle_test_brain_tailTwo(1,j,i) = all_median(1,j,2,i) - all_median(2,j,2,i) >= 0;
        shuffle_test_brain_tailOne(2,j,i) = all_median(1,j,2,i) - all_median(3,j,2,i) <= 0;
        shuffle_test_brain_tailTwo(2,j,i) = all_median(1,j,2,i) - all_median(3,j,2,i) >= 0;
        shuffle_test_brain_tailOne(3,j,i) = all_median(2,j,2,i) - all_median(3,j,2,i) <= 0;
        shuffle_test_brain_tailTwo(3,j,i) = all_median(2,j,2,i) - all_median(3,j,2,i) >= 0;
    end    
end

% get the median correction magnitudes
for i = 1:(size(summary_data,1)/2)
    summary_data{i,1} = cellfun(@median, summary_data{i,2}) - cellfun(@median, summary_data{i,1});
    summary_data{i,3} = cellfun(@median, summary_data{i,3}) - cellfun(@median, summary_data{i,2});
end
summary_data=summary_data(1:i,[1 3]);

% get the raw p-values of stats tests
shuffle_test_laser_tailOne = mean(shuffle_test_laser_tailOne,3);
shuffle_test_laser_tailTwo = mean(shuffle_test_laser_tailTwo,3);
shuffle_test_brain_tailOne = mean(shuffle_test_brain_tailOne,3);
shuffle_test_brain_tailTwo = mean(shuffle_test_brain_tailTwo,3);

plot_data_indiv = [];
plot_data_summary = [];
figure('Name', fig_name,'NumberTitle','off'); 
hold on
laser_interval = 0.5*[-1,1];
hemp_interval = 0.5;
angle_sig_laser = nan(1,size(summary_data,2)*size(summary_data,1));
angle_sig_contraSC = nan(1,size(summary_data,2)*size(summary_data,1));
stats_reporting={};
for i = 1:size(summary_data,2)
    for j = 1:size(summary_data,1)
        for k = 1:size(summary_data{j,i},2)
            plot_data = summary_data{j,i}(:,k);
            scatter((((i-1)*(size(summary_data,1)+hemp_interval)+j)*2+laser_interval(k))*ones(1,numel(plot_data)),plot_data, [], [.7 .7 .7],"filled")
            errorbar((((i-1)*(size(summary_data,1)+hemp_interval)+j)*2+laser_interval(k)),median(plot_data), median(plot_data)-prctile(all_median(j,i,k,:),25),prctile(all_median(j,i,k,:),75)-median(plot_data),'Marker','.','MarkerSize',40)
            stats_reporting=horzcat(stats_reporting,[median(plot_data),prctile(all_median(j,i,k,:),25),prctile(all_median(j,i,k,:),75)]);
            plot_data_summary = [plot_data_summary;[median(plot_data),prctile(all_median(j,i,k,:),25),prctile(all_median(j,i,k,:),75)]];
            plot_data_indiv = [(((i-1)*(size(summary_data,1)+hemp_interval)+j)*2+laser_interval(k)),plot_data'];
        end
        % if you need rank-sum tests
        [angle_sig_laser((i-1)*size(summary_data,1)+j), ~, ~] = ranksum(summary_data{j,i}(:,1),summary_data{j,i}(:,2),'alpha', 0.05, 'tail', 'both');
        [angle_sig_contraSC((i-1)*size(summary_data,1)+j), ~, ~] = ranksum(summary_data{j,i}(:,2),summary_data{j,2}(:,2),'alpha', 0.05, 'tail', 'both');
    end
end
xlim([0 15])
ylim(ylim_val)
end