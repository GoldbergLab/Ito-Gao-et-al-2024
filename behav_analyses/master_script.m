%% Analysis code for producing figures in Ito, Gau et. al (2024)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Enter the path to the root data directory (the one that contains all the .mat data files)
data_root = 'Z:\bsi8\Ito, Gao et. al (2024)\Revision_2\Open data repo\datasets';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For each section below, adjust the parameters to match each desired
% figure panel, and run this script to analyze the data and generate the
% figure.

%% Fig 1d-f, Extd Fig 2g-h, Fig 5b: Analysis of bilateral optogenetic inhibition or behavior only (L2 or L4 recentering)
fig_name = "Fig 1d-f, Extd Fig 2g-h, Fig 5b: Analysis of bilateral optogenetic inhibition or behavior only (L2 or L4 recentering)";

% Parameter options:
%
%   figure panel            data_set                    lick_nums   recentering    ylabel_input     notes
%   ------------------------------------------------------------------------------------------------------------------
%   Fig 1d left             behavior_2D_data.mat        1           0              contactSite      example mouse is #9 in the dataset
%   Fig 1d center           "                           2           0              "                "
%   Fig 1d right            "                           3           0              "                "
%   Fig 1e left             behavior_2D_data.mat        1           0              protAngle        example mouse is #1 in the dataset
%   Fig 1e center           behavior_2D_data.mat        2           0              "                "
%   Fig 1e right            behavior_2D_data.mat        3           0              "                "
%   Fig 1f left             ALM_bilateral_data.mat      3           0              protAngle
%   Fig 1f center           TJM1_bilateral_data.mat     3           0              "
%   Fig 1f right            TJS1_bilateral_data.mat     3           0              "
%   Fig 4b top left         behavior_recenter_data.mat  3           1              contactSite      example mouse is #2 in the dataset
%   Fig 4b top center       behavior_recenter_data.mat  4           1              "                "
%   Fig 4b top right        behavior_recenter_data.mat  5           1              "                "
%   Fig 4b bottom left      behavior_recenter_data.mat  3           1              protAngle        "
%   Fig 4b bottom center    behavior_recenter_data.mat  4           1              "                "
%   Fig 4b bottom right     behavior_recenter_data.mat  5           1              "                "
%   Extd Fig 2g left        ALM_bilateral_data.mat      3           0              duration
%   Extd Fig 2g center      TJM1_bilateral_data.mat     3           0              "
%   Extd Fig 2g right       TJS1_bilateral_data.mat     3           0              "
%   Extd Fig 2h left        ALM_bilateral_data.mat      3           0              pathlength
%   Extd Fig 2h center      TJM1_bilateral_data.mat     3           0              "
%   Extd Fig 2h right       TJS1_bilateral_data.mat     3           0              "

% Which data set to use (see table above)
data_set = 'behavior_2D_data.mat';
% Which lick number to analyze (see table above)
lick_num = 1;
% Non-recentering (0) or recentering (1)?
recentering = 0;
% Which kinematic variable to analyze?
ylabel_input = 'protAngle';

% modify the following variables to adjust figure paremeters
ylim_val = [-20 20]; % adjust plot y-axis of summary plots, change based on the variable analyzed
min_max_angle = [-30 30]; % range of histogram, change based on the variable analyzed
plot_individual = 0; % whether to plot for each mouse


% leave the following to default
contact_dur_opto = 0; % whether laser is only on during contact, deafult as 0 for analyses used in the paper
lick_index_contact_num = 1; % which lick index system to use, set to 1 by default
L2contact_only = 1; % whether to exclude licks without L2 contacts, set to 1 by default
mm_pix = 0.06; % pixel to mm conversion
bin_size = 2; 

% Construct full path
data_set = fullfile(data_root, data_set);

[plot_data_summary,plot_data_indiv,shuffle_test_laser_tailOne,shuffle_test_laser_tailTwo,shuffle_test_position_tailOne,shuffle_test_position_tailTwo] = ...
    plot_bilateral_opto(data_set,lick_num,recentering,ylabel_input,ylim_val,min_max_angle,plot_individual,contact_dur_opto,lick_index_contact_num,L2contact_only,mm_pix,bin_size, fig_name);

% for shuffle_test_laser, each p-value is for a comparison between laser off vs
% on, for left, center, and rights trials respectively
% for shuffle_test_position, first row is for laser off, and second row is
% for laser on. In each row, each p-value is for a comparison against the
% center condition, for left, center, and rights trials respectively

% for example histogram in Fig 1d, run the with
% 'behavior_2D_data' and ylabel_input = 'contactSite' and pick #9 (L1-3).
% for example histogram in Fig 1e, run the with
% 'behavior_2D_data' and ylabel_input = 'protAngle' and pick #7 (L1-3).

%% Fig 2b, Extd Fig 2b-d: compare protrusion probability between laser-on and laser-off
fig_name = "Fig 2b, Extd Fig 2b-d: compare protrusion probability between laser-on and laser-off";

% Parameter options:
%
%   figure panel    data_set                    
%   ------------------------------------------------------------------------------------------------------------------
%   Extd Fig 2b     ALM_bilateral_data.mat
%   Extd Fig 2c     TJM1_bilateral_data.mat
%   Extd Fig 2d     TJS1_bilateral_data.mat
%   Fig 2b          SC_bilateral_data.mat 

% load one dataset to analyze:
data_set = 'SC_bilateral_data.mat';


% leave the following to default (will analyze all three conditions (left/center/right))
trial_type = [1,2,3];
% 1 = protrusion probability, 2 = contact probability
probability_type = 1;

% Construct full path
data_set = fullfile(data_root, data_set);

[plot_data_summary,plot_data_indiv,shuffle_test_laser_tailOne,shuffle_test_laser_tailTwo] = plot_prot_probability(data_set,trial_type,probability_type, fig_name);

% for shuffle_test_laser, each p-value is for a comparison between laser off vs
% on, for left, center, and rights trials respectively

%% Fig 2e, Extd Fig 5c-f, Analysis of unilateral optogenetic inhibition by brain regions (combined left&right manipulations)
fig_name = "Fig 2e, Extd Fig 5c-f, Analysis of unilateral optogenetic inhibition by brain regions (combined left&right manipulations)";

% Parameter options:
%
%   figure panel    data_set_all                                        ylabel_input
%   ------------------------------------------------------------------------------------------------------------------
%   Fig 1e left         {'ALM_left_data.mat','ALM_right_data.mat'}      protAngle
%   Fig 1e center       {'TJM1_left_data.mat','TJM1_right_data.mat'}    "
%   Fig 1e right        {'SC_left_data.mat','SC_right_data.mat'}        "
%   Extd Fig 5c left    {'ALM_left_data.mat','ALM_right_data.mat'}      duration
%   Extd Fig 5c left    {'TJM1_left_data.mat','TJM1_right_data.mat'}    "
%   Extd Fig 5c left    {'SC_left_data.mat','SC_right_data.mat'}        "
%   Extd Fig 5d left    {'ALM_left_data.mat','ALM_right_data.mat'}      pathlength
%   Extd Fig 5d left    {'TJM1_left_data.mat','TJM1_right_data.mat'}    "
%   Extd Fig 5d left    {'SC_left_data.mat','SC_right_data.mat'}        "
%   Extd Fig 5e left    {'ALM_left_data.mat','ALM_right_data.mat'}      maxSpeed
%   Extd Fig 5e left    {'TJM1_left_data.mat','TJM1_right_data.mat'}    "
%   Extd Fig 5e left    {'SC_left_data.mat','SC_right_data.mat'}        "
%   Extd Fig 5f left    {'ALM_left_data.mat','ALM_right_data.mat'}      accPeakNum
%   Extd Fig 5f left    {'TJM1_left_data.mat','TJM1_right_data.mat'}    "
%   Extd Fig 5f left    {'SC_left_data.mat','SC_right_data.mat'}        "

% load one pair of datasets to analyze:
data_set_all = {'SC_vgat_CHR2_left_data.mat', 'SC_vgat_CHR2_right_data.mat'}; 
% select which kinematic variable to analyze, choose from: protAngle, latDisp, duration, pathlength, maxSpeed, accPeakNum, InterLickInterval, contactON_lickOFF, contactON_protOFF, contactSite, retractAngle:
ylabel_input = 'protAngle'; 

% modify the following variables to adjust figure paremeters
% adjust plot y-axis of summary plots, change based on the variable analyzed
ylim_val=[-20 20]; 


% leave the following to default
lick_num = 3; % L2: where the nick happens, L3 where lick re-aiming happens, L4 where recentering nick happens, L5: where recentering happens, select only one lick here, if you want to compare across licks, use the code in the last section
recentering = 0; % whether this is to compare recenter vs non-recenter trials (1=yes, 2=no)
lick_index_contact_num = 1; % which lick index system to use, set to 1 by default
L2contact_only = 1; % whether to exclude licks without L2 contacts, set to 1 by default
mm_pix = 0.06; % pixel to mm conversion

% Construct full paths
data_set_all = cellfun(@(data_set)fullfile(data_root, data_set), data_set_all, 'UniformOutput', false);

angle_sum = combined_opto_data(data_set_all,lick_num,recentering,ylabel_input,lick_index_contact_num,L2contact_only,mm_pix);
[plot_data_summary,plot_data_indiv,shuffle_test_laser_tailOne,shuffle_test_laser_tailTwo] = plot_unilateral_opto(angle_sum,ylabel_input,ylim_val, fig_name);

% for shuffle_test_laser, each p-value is for a comparison between laser off vs
% on, for left, center, and rights trials respectively

% for example histogram in Fig 2c, run the previous section with
% 'TJM1_right_data.mat' and pick #3 for left column, run with 'TJM1_left_data.mat' 
% and pick #3 for center column, run with 'SC_left_data.mat' 
% and pick #4 for center column
%% Fig 2f: calculate correction magnitude of brain region x sides
fig_name = "Fig 2f: calculate correction magnitude of brain region x sides";

% Parameter options:
%
% None

% combine datasets to analyze:
data_set_all = {'ALM_left_data.mat', 'TJM1_left_data.mat',...
    'SC_left_data.mat','ALM_right_data.mat'...
    'TJM1_right_data.mat','SC_right_data.mat'}; 

% modify the following variables to adjust figure paremeters
ylim_val=[0 15]; % adjust plot y-axis of summary plots, change based on the variable analyzed


% leave the following to default
% which kinematic variable to analyze, choose from: protAngle, latDisp, duration, pathlength, maxSpeed, accPeakNum, InterLickInterval, contactON_lickOFF, contactON_protOFF, contactSite, retractAngle
ylabel_input = 'protAngle';
% L2: where the nick happens, L3 where lick re-aiming happens, L4 where recentering nick happens, L5: where recentering happens, select only one lick here, if you want to compare across licks, use the code in the last section
lick_num = 3; 
% whether this is to compare recenter vs non-recenter trials (1=yes, 2=no)
recentering = 0; 
% which lick index system to use, set to 1 by default
lick_index_contact_num = 1; 
% whether to exclude licks without L2 contacts, set to 1 by default
L2contact_only = 1; 
% pixel to mm conversion
mm_pix = 0.06; 

% Construct full paths
data_set_all = cellfun(@(data_set)fullfile(data_root, data_set), data_set_all, 'UniformOutput', false);

angle_sum = combined_opto_data(data_set_all,lick_num,recentering,ylabel_input,lick_index_contact_num,L2contact_only,mm_pix);
[plot_data_summary,plot_data_indiv,shuffle_test_laser_tailOne,shuffle_test_laser_tailTwo,shuffle_test_brain_tailOne,shuffle_test_brain_tailTwo] = plot_correction_magnitude(angle_sum,ylim_val, fig_name);
% for shuffle_test_laser, first column is for ipsilateral corrections and 
% the second column contralateral corrections. First row is ALM, second TJM1
% and third SC. Eachp-value is for a comparison between laser off vs
% on.
% for shuffle_test_brain, first column is for ipsilateral corrections and 
% the second column contralateral corrections. First row is ALM vs TJM1, 
% second ALM vs SC, and third TJM1 vs SC

%% Fig 5b: cumulative protrusion probability with stimulation
fig_name = "Fig 5b: cumulative protrusion probability with stimulation";

% Parameter options:
%
% None

% Notes:
% In each mouse's dataset for stimulation (t_stats_combined(xx).data/video_descriptor), rows 1-5 represent stimulation sessions
% with AP of -2.9(not used), -3,1, -3.4. -3.7, -3.9, respectively. Columns
% 1-3 represent ML of +-1.3,+-1.1,+-0.9(not used), respectively for AP
% -2.9-(-3.1); ML of +-1.6,+-1.3,+-1.1(not used), respectively for AP
% -3.3-(-3.9). Video descriptor is used to figure the number of trials in
% each session and which of them had laser on

% load one dataset to analyze:
data_set = 'stim_uncued_data.mat';
% which trial type to analyze, select from 'ctrl','cued','uncued'
stim_type = 'uncued';
animal_ID = [1,3,5,6,7,8]; % for this analysis, we only considered sessions with activations of the same 400ms durations


% leave the following to default:
start_column=1; % most lateral site to analyze. Based on histology, all AP locations have at most two different sites in the ML direction. They are always input to the first two columns.
end_column=2; % most medial site to analyze
start_row=2; % most anterior site to analyze. In this code, we started from -3.1AP as -2.9AP sites were outsite SC.
end_row=5; % most analyze site to combine

% Construct full path
data_set = fullfile(data_root, data_set);

[plot_summary,plot_indiv] = plot_stimulation_xsite_cumulProb(data_set,stim_type,animal_ID,start_column,end_column,start_row,end_row, fig_name);

%% Fig 5c, Extd Fig 10d: example polar plots showing protrusion endpoints
fig_name = "Fig 5c, Extd Fig 10d: example polar plots showing protrusion endpoints";

% Parameter options:
%
%   figure panel        data_set                animal_num
%   ------------------------------------------------------------------------------------------------------------------
%   Fig 5c left         stim_uncued_data.mat    4
%   Fig 5c right        "                       7
%   Extd Fig 10d left   stim_cued_data.mat      5
%   Extd Fig 10d right  "                       8

% load one dataset to analyze:
data_set = 'stim_cued_data.mat';
% Specify which animal to analyze
animal_num = 8;
% Note that for left examples, plot must be horizontally flipped


% 1=left, -1=right, this array is for all 8 mice used here
activationSide=[-1,1,1,1,1,-1,-1,-1]; 

% Construct full path
data_set = fullfile(data_root, data_set);

plot_stimulation_polar_plot(data_set,animal_num,activationSide, fig_name)

%% Fig 5d&f, Extd Fig 10e&g: compare across stimulation sites (except for comparing protrusion probability)
fig_name = "Fig 5d&f, Extd Fig 10e&g: compare across stimulation sites (except for comparing protrusion probability)";

% Parameter options:
%
%   figure panel        data_set                ylabel_input    laser_type
%   ------------------------------------------------------------------------------------------------------------------
%   Fig 5d              stim_uncued_data.mat    protAngle       1
%   Fig 5f              "                       protLatency     1
%   Extd Fig 10e        stim_cued_data.mat      protAngle       1
%   Extd Fig 10g left   "                       protLatency     1
%   Extd Fig 10g right  "                       protLatency     0

% load one dataset to analyze:
data_set = 'stim_uncued_data.mat';
% which kinematic variable to analyze, choose from: protAngle, protLatency
ylabel_input = 'protLatency'; 
laser_type = 1; % 1=with laser, 0=no laser


% leave the following to default:
activationSide=[-1,1,1,1,1,-1,-1,-1];   %1=left, -1=right, this array is for all 8 mice used here
lick_num=1;
site_AP=[2.9 3.1,3.3,3.7,3.9];          % all AP coordinates with data collected, note -2.9 isn't used in the analysis because it may be outside SC
start_column=1;                         % most lateral site to analyze. Based on histology, all AP locations have at most two different sites in the ML direction. They are always input to the first two columns.
end_column=2;                           % most medial site to analyze
start_row=2;                            % most anterior site to analyze. In this code, we started from -3.1AP as -2.9AP sites were outsite SC.
end_row=5;                              % most analyze site to combine

% Construct full path
data_set = fullfile(data_root, data_set);

[plot_data_summary,plot_data_indiv,p_kw,c_kw,corr_p,shuffle_test_p,observed_r_sqr,shuffled_r_sqr] = plot_stimulation_xsite_comparison(data_set,ylabel_input,activationSide,lick_num,laser_type,site_AP,start_column,end_column,start_row,end_row, fig_name);

% p_kw and c_kw contains p-values of the Kruskal-Wallis test and its
% post-hoc test respectively

%% Fig 5e, Extd Fig 10f: compare protrusion probability across stimulation sites
fig_name = "Fig 5e, Extd Fig 10f: compare protrusion probability across stimulation sites";

% Parameter options:
%
%   figure panel        data_set                stim_type ylabel_input    laser_type
%   ------------------------------------------------------------------------------------------------------------------
%   Fig 5e              stim_uncued_data.mat    uncued
%   Extd Fig 10f left   stim_cued_data.mat      cued
%   Extd Fig 10f right  stim_cued_data.mat      ctrl

% load one dataset to analyze:
data_set = 'stim_cued_data.mat';
% which trial type to analyze, select from 'ctrl','cued','uncued'
stim_type = 'ctrl'; 


% leave the following to default:
start_column=1; % most lateral site to analyze. Based on histology, all AP locations have at most two different sites in the ML direction. They are always input to the first two columns.
end_column=2; % most medial site to analyze
start_row=2; % most anterior site to analyze. In this code, we started from -3.1AP as -2.9AP sites were outsite SC.
end_row=5; % most analyze site to combine

% Construct full path
data_set = fullfile(data_root, data_set);

[plot_data_summary,plot_data_indiv,p_kw,c_kw] = plot_stimulation_xsite_protProb(data_set,stim_type,start_column,end_column,start_row,end_row, fig_name);
% p_kw and c_kw contains p-values of the Kruskal-Wallis test and its
% post-hoc test respectively

%% Extd Fig 1, Extd Fig 4: compare two sets of sessions
fig_name = "Extd Fig 1, Extd Fig 4: compare two sets of sessions";

% Parameter options:
%
%   figure panel        data_set_one            data_set_two            lick_num    true_L1     ylabel_input    ylim_val    notes
%   --------------------------------------------------------------------------------------------------------------------------------------------------------------
%   Extd Fig 1d left    L2_dispense_data.mat    L3_dispense_data.mat    1           0           protAngle       [-20 20]    example mouse is #7 in the dataset for Extd Fig 1c 
%   Extd Fig 1d cemter  L2_dispense_data.mat    L3_dispense_data.mat    2           0           "               "
%   Extd Fig 1d right   L2_dispense_data.mat    L3_dispense_data.mat    3           0           "               "
%   Extd Fig 4b         pre_lesion_data.mat     post_lesion_data.mat    3           0           protAngle       [-20 20]
%   Extd Fig 4d         "                       "                       1           1           reaction_time   [0 500]
%   Extd Fig 4e         "                       "                       1           1           duration        [0 100]
%   Extd Fig 4f         "                       "                       1           1           pathlength      [0 16]
%   Extd Fig 4g         "                       "                       1           1           maxSpeed        [0 400]
%   Extd Fig 4h         "                       "                       1           1           accPeakNum      [0 10]

% load two datasets to analyze:
data_set_one = 'pre_lesion_data.mat';
data_set_two = 'post_lesion_data.mat'; 

% select what to analyze:
lick_num = 1; % L2: where the nick happens, L3 where lick re-aiming happens, L4 where recentering nick happens, L5: where recentering happens, select only one lick here, if you want to compare across licks, use the code in the last section
true_L1 = 1; % 1: the first lick after the cue; 0: first lick that makes spout contact
ylabel_input = 'accPeakNum'; % which kinematic variable to analyze, choose from: protAngle, latDisp, duration, pathlength, maxSpeed, accPeakNum, InterLickInterval, contactON_lickOFF, contactON_protOFF, contactSite, retractAngle, reaction_time

% modify the following variables to adjust figure paremeters
ylim_val=[0 10]; % adjust plot y-axis of summary plots, change based on the variable analyzed
min_max_angle = [-30 30]; % range of histogram, change based on the variable analyzed
plot_individual = 0; % whether to plot for each mouse


% leave the following to default
contact_dur_opto = 0; % whether laser is only on during contact, deafult as 0 for analyses used in the paper
lick_index_contact_num = 1; % which lick index system to use, set to 1 by default
L2contact_only = 1; % whether to exclude licks without L2 contacts, set to 1 by default
mm_pix = 0.06; % pixel to mm conversion
bin_size = 2; % bin size of histogram

% Construct full paths
data_set_one = fullfile(data_root, data_set_one);
data_set_two = fullfile(data_root, data_set_two);

[plot_data_summary,plot_data_indiv,shuffle_test_lesion_tailOne,shuffle_test_lesion_tailTwo,shuffle_test_position_tailOne,shuffle_test_position_tailTwo,shuffle_test_lesion_combined_tailOne,shuffle_test_lesion_combined_tailTwo,plot_data_summary_combined,plot_data_indiv_combined] = ...
    plot_xsession_comparison(data_set_one,data_set_two,lick_num,true_L1,ylabel_input,ylim_val,min_max_angle,plot_individual,lick_index_contact_num,L2contact_only,mm_pix,bin_size, fig_name);

% for shuffle_test_lesion, each p-value is for a comparison between pre- vs
% post-lesion, for left, center, and rights trials respectively
% for shuffle_test_position, first row is for pre-lesion, and second row is
% for post-lesion. In each row, each p-value is for a comparison against the
% center condition, for left, center, and rights trials respectively
% shuffle_test_lesion_combined is for a pre- vs post-lesion comparison
% with data pooled across all trial types

%% Extd Fig 2a: compare L1>L2 vs L2>L3 interlick intervals
fig_name = "Extd Fig 2a: compare L1>L2 vs L2>L3 interlick intervals";

% Parameter options:
%
% None

% load one dataset to analyze:
data_set = 'behavior_2D_data.mat'; 
% select what to analyze:
ylabel_input = 'InterLickInterval'; % which kinematic variable to analyze, choose from: protAngle, latDisp, duration, pathlength, maxSpeed, accPeakNum, InterLickInterval, contactON_lickOFF, contactON_protOFF, contactSite, retractAngle
lick_num = [2,3]; % select which two licks to compare: L2: where the nicks happens, L3 where lick re-aiming happens, L4 where recentering nick happens, L5: where recentering happens
laser_on = 0; % this code only allows you to look at laser OR no-laser trials
recentering = 0; % whether to compare recenter vs non-recenter trials (1=yes, 2=no)
% modify the following variables to adjust figure paremeters
ylim_val=[0 200]; % adjust plot y-axis of summary plots based on the variable analyzed

% leave the following to default
plot_individual = 0; % whether to plot for each mouse
lick_index_contact_num = 1; % which lick index system to use, set to 1 by default
L2contact_only = 1; % whether to exclude licks without L2 contacts, set to 1 by default
mm_pix = 0.06; % pixel to mm conversion
min_max_angle = [-30 30]; % range of histogram
bin_size = 2; % bin size of histogram

% Construct full path
data_set = fullfile(data_root, data_set);

[plot_data_summary,plot_data_indiv,shuffle_test_lick_tailOne,shuffle_test_lick_tailTwo,shuffle_test_position_tailOne,shuffle_test_position_tailTwo] = plot_xlick_comparison(data_set,lick_num,laser_on,recentering,ylabel_input,ylim_val,min_max_angle,plot_individual,contact_dur_opto,lick_index_contact_num,L2contact_only,mm_pix,bin_size, fig_name);

% for shuffle_test_lick, each p-value is for a comparison between the two
% licks, for left, center, and rights trials respectively
% for shuffle_test_position, first row is for the first lick chosen, and 
% second row is for the second lick chosen. In each row, each p-value is 
% for a comparison against the center condition, for left, center, and 
% rights trials respectively

%% Extd Fig 2e-f, Extd Fig 5a-b: compare protrusion probability between laser-on and laser-off, for each trial type
fig_name = "Extd Fig 2e-f, Extd Fig 5a-b: compare protrusion probability between laser-on and laser-off, for each trial type";

% Parameter options:
%
%   figure panel        data_set                                        probability_type
%   ------------------------------------------------------------------------------------------------------------------
%   Extd Fig 2e left    ALM_bilateral_data.mat                          1
%   Extd Fig 2e center  TJM1_bilateral_data.mat                         1
%   Extd Fig 2e right   TJS1_bilateral_data.mat                         1
%   Extd Fig 2f left    ALM_bilateral_data.mat                          2
%   Extd Fig 2f center  TJM1_bilateral_data.mat                         2
%   Extd Fig 2f right   TJS1_bilateral_data.mat                         2
%   Extd Fig 5a left    {'ALM_left_data.mat','ALM_right_data.mat'}      1
%   Extd Fig 5a center  {'TJM1_left_data.mat','TJM1_right_data.mat'}    1
%   Extd Fig 5a right   {'SC_left_data.mat','SC_right_data.mat'}        1
%   Extd Fig 5b left    {'ALM_left_data.mat','ALM_right_data.mat'}      2
%   Extd Fig 5b center  {'TJM1_left_data.mat','TJM1_right_data.mat'}    2
%   Extd Fig 5b right   {'SC_left_data.mat','SC_right_data.mat'}        2

% load one dataset to analyze:
data_set = {'ALM_left_data.mat','ALM_right_data.mat'};


% leave the following to default
trial_type = [1,2,3]; % default setting will analyze all three conditions (left/center/right)
probability_type = 2; % 1 = contact probability, 2 = CSM probability

% Construct full paths
if iscell(data_set)
    data_set = cellfun(@(data_set)fullfile(data_root, data_set), data_set, 'UniformOutput', false);
else
    data_set = fullfile(data_root, data_set);
end

[plot_data_summary,plot_data_indiv,shuffle_test_laser_tailOne,shuffle_test_laser_tailTwo] = plot_probability_trial_type(data_set,trial_type,probability_type, fig_name);

% for shuffle_test_laser, each p-value is for a comparison between laser off vs
% on, for left, center, and rights trials respectively

%% Extd Fig 4c: compare two sets of sessions for trial initiation
fig_name = "Extd Fig 4c: compare two sets of sessions for trial initiation";

% Parameter options:
%
% None

% Extd Fig 4c, data_set_one = 'pre_lesion_data.mat'; data_set_two = 'post_lesion_data.mat';
data_set_one = 'pre_lesion_data.mat';
data_set_two = 'post_lesion_data.mat'; 

% Construct full paths
data_set_one = fullfile(data_root, data_set_one);
data_set_two = fullfile(data_root, data_set_two);

[shuffle_test_lesion_combined_tailOne,shuffle_test_lesion_combined_tailTwo,plot_data_summary_combined,plot_data_indiv_combined] = plot_xsession_comparison_initiation(data_set_one,data_set_two, fig_name);

%% Extd Fig 7b-f example correlation of contact angle vs protrusion angle change
fig_name = "Extd Fig 7b-f example correlation of contact angle vs protrusion angle change";

% Parameter options:
%
%   figure panel    data_set                    lick_1  lick_2  indiv_trial_type
%   ------------------------------------------------------------------------------------------------------------------
%   Extd Fig 7b     example_correlation.mat     1       2       [1, 2, 3]
%   Extd Fig 7c     "                           2       3       [1, 2, 3]
%   Extd Fig 7d     "                           2       3       1
%   Extd Fig 7e     "                           2       3       2
%   Extd Fig 7f     "                           2       3       3

% load one file to analyze:
data_set = 'example_correlation.mat';
% define licks numbers to analyze
lick_1 = 2; 
lick_2 = 3;
% which conditions to plot individually
indiv_trial_type = [1,2,3]; 


% Construct full path
data_set = fullfile(data_root, data_set);

plot_angle_correlation_example(data_set,lick_1,lick_2,indiv_trial_type, fig_name)

%% Extd Fig 7b-f summary correlation of contact angle vs protrusion angle change
fig_name = "Extd Fig 7b-f summary correlation of contact angle vs protrusion angle change";

% Parameter options:
%
%   figure panel    data_set                    lick_1  lick_2  indiv_trial_type
%   ------------------------------------------------------------------------------------------------------------------
%   Extd Fig 7b     behavior_correlation_data.mat  1       2       [1, 2, 3]
%   Extd Fig 7c     behavior_correlation_data.mat  2       3       [1, 2, 3]
%   Extd Fig 7d     behavior_correlation_data.mat  2       3       1
%   Extd Fig 7e     behavior_correlation_data.mat  2       3       2
%   Extd Fig 7f     behavior_correlation_data.mat  2       3       3

% load one dataset to analyze:
data_set = 'behavior_correlation_data.mat';
% define licks to analyze
lick_1 = 1; 
lick_2 = 2;
% which conditions to plot individually
indiv_trial_type = [1,2,3]; 

% Construct full path
data_set = fullfile(data_root, data_set);

[plot_data_summary,plot_data_indiv,p_value_tailOne_summary,p_value_tailTwo_summary] = plot_angle_correlation_summary(data_set,lick_1,lick_2,indiv_trial_type, fig_name); %#ok<*ASGLU> 

%% Extd Fig 10h-m: kinematic comparisons across cued stimulation, uncued stimulation, and cue only conditions
fig_name = "Extd Fig 10h-m: kinematic comparisons across cued stimulation, uncued stimulation, and cue only conditions";

% Parameter options:
%
%   figure panel    data_set_cued           data_set_uncued         ylabel_input
%   ------------------------------------------------------------------------------------------------------------------
%   Extd Fig 10h    stim_cued_data.mat    stim_uncued_data.mat    protLatency
%   Extd Fig 10i    stim_cued_data.mat    stim_uncued_data.mat    duration
%   Extd Fig 10j    stim_cued_data.mat    stim_uncued_data.mat    InterLickInterval
%   Extd Fig 10k    stim_cued_data.mat    stim_uncued_data.mat    pathlength
%   Extd Fig 10l    stim_cued_data.mat    stim_uncued_data.mat    maxSpeed
%   Extd Fig 10m    stim_cued_data.mat    stim_uncued_data.mat    accPeakNum

% load two dataseta to analyze:
data_set_cued = 'stim_cued_data.mat';
data_set_uncued = 'stim_uncued_data.mat';

% which kinematic variable to analyze, choose from: protLatency, duration, InterLickInterval, pathlength, maxSpeed, accPeakNum
ylabel_input = 'protLatency'; 

% modify the following variables to adjust figure paremeters
% adjust plot y-axis of summary plots, change based on the variable analyzed
ylim_val=[0 250]; 


% leave the following to default:
lick_num=1;

% Construct full paths
data_set_cued = fullfile(data_root, data_set_cued);
data_set_uncued = fullfile(data_root, data_set_uncued);

[plot_data_summary,plot_data_indiv,shuffled_test_results] = plot_stimulation_kinematics_comparison(data_set_cued,data_set_uncued,ylabel_input,lick_num,ylim_val, fig_name);
% for shuffle_test_results, first column is cue only vs stim. w/ cue,
% second column is cue only vs stim. only, third column is stim. w/ cue vs
% stim. only

