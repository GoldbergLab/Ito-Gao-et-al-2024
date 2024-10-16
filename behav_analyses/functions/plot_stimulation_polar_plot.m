function plot_stimulation_polar_plot(data_set,animal_num,activationSide, fig_name)

load(data_set);
data_indiv = t_stats_combined(animal_num).data;
activationSide = activationSide(animal_num);

% the following should be kept default
theta_offset=+pi/2;
colors={'yellow','red','green','magenta'};
start_row=2; %first row to use
end_row = 5; %last row to use
mm_pix = 0.06; % pixel to mm conversion
lick_num=1;
laser_type=1; % 1=with laser, 0=no laser
data_x=cell(1,end_row);
data_y=cell(1,end_row);

for i=start_row:end_row
    temp_x = [];
    temp_y = [];
    % for example polar plots, we are using only sites with ML = 1.6
    if i == 2 % at AP = -3.1, the first column is ML1.6
        col_num = 1;
    elseif i > 2 % at AP more posterior than -3.1, the second column is ML1.6
        col_num = 2;
    end
    if ~isempty(data_indiv{i,col_num})
        t_stats = data_indiv{i,col_num};
        for p=1:numel(t_stats)
            x = t_stats(p).tip_xf*activationSide;
            y = t_stats(p).tip_yf;
            t = min(find(x==max(x)));
            index = t_stats(p).prot_ind;
            if t_stats(p).lick_index == lick_num && t_stats(p).laser_trial == 1 &&t_stats(p).laser_2D == laser_type && t_stats(p).time_rel_cue>0
                temp_x = [temp_x, x(t)*mm_pix];
                temp_y = [temp_y, y(t)*mm_pix];
            end
        end
    end
    data_x{i}=temp_x;
    data_y{i}=temp_y;
end

figure('Name', fig_name,'NumberTitle','off'); 
loop_count=1;
for i=start_row:end_row
    theta=atan(data_x{i}./data_y{i});
    rho=sqrt(data_x{i}.^2+data_y{i}.^2);
    p=polarplot(theta+theta_offset,rho,"x");
    p.Color=colors{loop_count};
    hold on
    index=find(theta==median(theta));
    p=polarplot([0,median(theta)+theta_offset],[0,median(rho)],"--");
    p.Color=colors{loop_count};
    loop_count=loop_count+1;
end
rticks([0 4 8])
rlim([0 8])
thetalim([0 180])
thetaticks([0:45:180])
thetaticklabels({'90','45','0','-45','-90'})

end