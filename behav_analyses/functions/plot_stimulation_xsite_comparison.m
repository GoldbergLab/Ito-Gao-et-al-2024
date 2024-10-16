function [plot_data_summary,plot_data_indiv,p_kw,c_kw,corr_p,shuffle_test_p,observed_r_sqr,shuffled_r_sqr] = plot_stimulation_xsite_comparison(data_set,ylabel_input,activationSide,lick_num,laser_type,site_AP,start_column,end_column,start_row,end_row, fig_name)

load(data_set);
gridSize=[5,3]; % see description above
loop_num=1000; % for shuffles and bootstraps
data_raw=cell(gridSize);

for i=1:gridSize(1)
    for j=1:gridSize(2)
        temp_data=cell(1,numel(t_stats_combined));
        for k=1:numel(t_stats_combined)
            data_indiv = t_stats_combined(k).data;
            if ~isempty(data_indiv{i,j})
                t_stats = data_indiv{i,j};
                for p=1:numel(t_stats)
                    x = t_stats(p).tip_xf*activationSide(k);
                    y = t_stats(p).tip_yf;
                    t = min(find(x==max(x)));
                    index = t_stats(p).prot_ind;
                    if t_stats(p).lick_index == lick_num && t_stats(p).laser_trial == laser_type &&t_stats(p).laser_2D == laser_type && t_stats(p).time_rel_cue>0
                        if strcmp(ylabel_input,'protAngle')
                            temp_data{k}=[temp_data{k},atan((x(t))/y(t))*180/pi];
                        elseif strcmp(ylabel_input,'protLatency')
                            temp_data{k}=[temp_data{k},t_stats(p).time_rel_cue];
                        end
                    end
                end
            end
            data_raw{i,j}=temp_data;
        end
    end
end

% combine data with the same AP for each mouse
AP_data=cell(1,gridSize(1));
AP_indiv=cell(1,gridSize(1));
AP_median=NaN(1,gridSize(1));
for i=1:gridSize(1)
    temp_AP_data=cell(1,numel(t_stats_combined));
    for j=1:numel(t_stats_combined)
        temp_data=[];
        for k=start_column:end_column
            temp_data=[temp_data,data_raw{i,k}{j}];
        end
        temp_AP_data{j}=temp_data;
    end
    AP_data{i}=temp_AP_data;
    AP_indiv{i}=cellfun(@median,AP_data{i});
    AP_median(i)=median(AP_indiv{i},'omitnan');
end

%get bootstrapped IQR and correlation coefficients
AP_median_BS = nan(loop_num,numel(AP_data));
corr_BS = nan(loop_num,1);
for m=1:loop_num
    data_BS_xAnimal=datasample([1:numel(AP_data{1,1})],numel(AP_data{1,1}));
    median_BS_xAnimal=nan(1,numel(data_BS_xAnimal));
    stats_data = [];
    site_data = [];
    for n=1:numel(AP_data)
        data_temp = AP_data{n}(data_BS_xAnimal); % get data from bs animals for a stim site
        data_temp = cellfun(@(x) datasample(x,numel(x)),data_temp,'UniformOutput',false); % sample each animal's data
        data_temp = cellfun(@median, data_temp); % get median of each mouse
        AP_median_BS(m,n) = median(data_temp, 'omitnan');
        nanIND = isnan(data_temp);
        if n >= start_row
            stats_data=[stats_data;data_temp(~nanIND)'];
            site_data=[site_data;site_AP(n)*ones(numel(data_temp(~nanIND)),1)];
        end
    end
    [corr_BS(m),~] = corr(site_data,stats_data, 'Type', 'Spearman');
end

AP_iqr_h=prctile(AP_median_BS,75,1);
AP_iqr_l=prctile(AP_median_BS,25,1);
AP_indiv=AP_indiv(start_row:end_row);
AP_median=AP_median(start_row:end_row);
AP_iqr_l=AP_iqr_l(start_row:end_row);
AP_iqr_h=AP_iqr_h(start_row:end_row);

% kruskal wallis test across sites
stats_data=[];
for i=1:size(AP_indiv,2)
    stats_data=[stats_data,AP_indiv{1,i}'];
end
[p_kw,tbl,stats] = kruskalwallis(stats_data);
c_kw = multcompare(stats,'Display','off');

% correlation test across sites
site_AP=site_AP(start_row:end_row);
test_data=AP_indiv';
stats_data=[];
site_data=[];
for i=1:size(test_data)
    temp_data=test_data{i};
    nanIND=isnan(temp_data);
    stats_data=[stats_data;temp_data(~nanIND)'];
    site_data=[site_data;site_AP(i)*ones(numel(temp_data(~nanIND)),1)];
end
[actual_r,corr_p] = corr(site_data,stats_data, 'Type', 'Spearman');
actual_r_sqr = actual_r^2;

% plot data across sites, with lines connecting data points from the same mouse
figure('Name', fig_name,'NumberTitle','off'); 
randomColors=rand(numel(t_stats_combined),3);
grayColor=[.7 .7 .7];
blkColor=[0 0 0];
xval= [3.1 3.4 3.7 3.9];
for k=1:numel(t_stats_combined)
    for i=1:size(AP_indiv,2)-1
        if ~isnan(AP_indiv{1,i}(k)) && ~isnan(AP_indiv{1,i+1}(k))
            plot([xval(i),xval(i+1)],[AP_indiv{1,i}(k),AP_indiv{1,i+1}(k)],'color',grayColor,'MarkerEdgeColor',grayColor,'MarkerFaceColor',grayColor,'LineWidth',1)
            hold on
        end
    end
end
for i=1:size(AP_median,2)-1
    plot([xval(i),xval(i+1)],[AP_median(1,i),AP_median(1,i+1)],'color',blkColor,'MarkerEdgeColor',blkColor,'MarkerFaceColor',blkColor,'LineWidth',3)
    hold on
    errorbar(xval,AP_median,AP_median - AP_iqr_l,AP_iqr_h - AP_median,'.','MarkerSize',12) 
end
xlim([min(xval)-0.075,max(xval)+0.075])
xticks(xval)
xticklabels({'-3.1AP','-3.4AP','-3.7AP','-3.9AP'})
xlabel('AP coordinates')
ylabel(ylabel_input)
box off

% plot data across sites, as dots only
figure('Name', fig_name,'NumberTitle','off'); 
hold on
for i=1:size(AP_median,2)
    errorbar(xval(i),AP_median(1,i),AP_median(1,i)-AP_iqr_l(i),AP_iqr_h(i)-AP_median(i));
    scatter(ones(1,numel(AP_indiv{1,i}))*xval(i),AP_indiv{1,i},[], [0 0 0],"filled");
end
xlim([min(xval)-0.075,max(xval)+0.075])
xticks(xval)
xticklabels({'-3.1AP','-3.4AP','-3.7AP','-3.9AP'})
xlabel('AP coordinates')
ylabel(ylabel_input)
box off

for i = 1:numel(AP_median)
    plot_data_summary(i,1) = AP_median(i);
    plot_data_summary(i,2) = AP_iqr_l(i);
    plot_data_summary(i,3) = AP_iqr_h(i);
    plot_data_indiv(i,:) = AP_indiv{i};
end

% get correlation coefficients from shuffles
r_rand = nan(1,loop_num);
shuffle_test_p = nan(1,loop_num);
for m=1:loop_num
    rand_data = cell(numel(AP_data{1,1}),numel(AP_data));
    for n = 1:numel(AP_data{1,1})
        data_temp = cellfun(@(x) x(n), AP_data);
        numel_temp = cellfun(@(x) numel(x), data_temp);
        data_temp = horzcat(data_temp{:});
        rand_ind = randperm(numel(data_temp));
        for p = 1:numel(AP_data)
            if numel_temp(p) > 0
                if p == 1
                    rand_data{n,p} = data_temp(rand_ind([1:numel_temp(p)]));
                else
                    rand_data{n,p} = data_temp(rand_ind([1+sum(numel_temp(1:p-1)):sum(numel_temp(1:p))]));
                end
            end
        end
        
    end
    
    rand_data = rand_data(:,start_row:end_row);
    rand_data = cellfun(@median, rand_data);
    stats_data=[];
    site_data=[];
    for i=1:size(rand_data,2)
        temp_data=rand_data(:,i);
        nanIND=isnan(temp_data);
        stats_data=[stats_data;temp_data(~nanIND)];
        site_data=[site_data;site_AP(i)*ones(numel(temp_data(~nanIND)),1)];
    end
    [r_rand(m),~] = corr(site_data,stats_data, 'Type', 'Spearman');
    shuffle_test_p(m) = r_rand(m)^2 >= actual_r_sqr;
end
shuffle_test_p = sum(shuffle_test_p)/numel(shuffle_test_p);

% plot bootstrapped correlation coefficents with shuffled ones
figure('Name', fig_name,'NumberTitle','off'); 
errorbar(2,median(r_rand)^2,median(r_rand)^2-prctile(r_rand,25)^2,prctile(r_rand,75)^2-median(r_rand)^2,'.')
xlim([0.5 2.5])
xticks([1,2])
xticklabels({'Obs.','Shuff.'})
ylim([-0.1 1])
ylabel('R^2-value')
yticks([-1,0,1])
hold on
errorbar(1,actual_r_sqr, actual_r_sqr - prctile(corr_BS,25)^2, prctile(corr_BS,75)^2 - actual_r_sqr, '.')
box off
observed_r_sqr = [actual_r_sqr, prctile(corr_BS,25)^2, prctile(corr_BS,75)^2];
shuffled_r_sqr = [median(r_rand)^2, prctile(r_rand,25)^2, prctile(r_rand,75)^2];

end