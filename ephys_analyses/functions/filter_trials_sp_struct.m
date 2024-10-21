function filter_ind = filter_trials_sp_struct(sp_struct, neuron, lick_index, lick_exist_lick_num,...
    sp_contact_lick_num, sp_contact2_absent_lick_num, contact_centroid_lick_num)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns indices corresponding to trials in sp_struct that
% fulfill all criteria imposed by the experimenter.
%
% neuron is the neuron index in sp_struct of interest
%
% lick_index is the lick index to be used to filter trials
%
% lick_exist_lick_num is an array of lick numbers to check whether a lick
% exists or not. e.g., [1 2 3 4 5]
%
% sp_contact_lick_num is an array of lick numbers to check whether a
% lick made spout contact or not. e.g., [1 2 3 4 5]
%
% sp_contact2_absent_lick_num is an array of lick numbers to check 
% whether there is a double-tap on a lick or not. e.g., [1 2 3 4 5]
%
% contact_centroid_lick_num is an array of lick numbers to check 
% whether contact_centroid was found for a given lick or not.
%
% Author:  Brendan Ito
% Email:   bsi8@cornell.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% create filter for trials where licks exist
if ~isempty(lick_exist_lick_num)
    lick_exist = filter_lick_exist(lick_index, lick_exist_lick_num);
else
    lick_exist = ones(size(lick_index));
end

% create filter for trials where licks (lick_num) made spout contact
if ~isempty(sp_contact_lick_num)
    sp_contact = sp_struct(neuron).sp_contact;
    sp_contact_exist = filter_sp_contact(sp_contact, lick_index, sp_contact_lick_num);
else
    sp_contact_exist = ones(size(lick_index));
end
    
% create filter for trials where licks (lick_num) did not double-tap
if ~isempty(sp_contact2_absent_lick_num)
    sp_contact2_absent = sp_struct(neuron).sp_contact2_absent;
    sp_contact2_absent_exist = filter_sp_contact2_absent(sp_contact2_absent, lick_index, sp_contact2_absent_lick_num);
else
    sp_contact2_absent_exist = ones(size(lick_index));
end

% create filter for trials where licks (lick_num) made spout contact
if ~isempty(contact_centroid_lick_num)
    contact_centroid = sp_struct(neuron).contact_centroid;
    contact_centroid_exist = filter_contact_centroid(contact_centroid, lick_index, contact_centroid_lick_num);
else
    contact_centroid_exist = ones(size(lick_index));
end

% combine all filters
filter_ind = lick_exist & sp_contact_exist & sp_contact2_absent_exist & contact_centroid_exist;
end

%%

function output = filter_sp_contact(sp_contact, lick_index, lick_nums)

% find trials where licks made contact with the spout
output = ones(1, numel(lick_index));
for j = 1:numel(lick_nums)
    lick_contact_exist_temp = cellfun(@(x, y) ~isnan(y(x == lick_nums(j))), lick_index, sp_contact, 'UniformOutput', false);
    empty_index = cellfun('isempty', lick_contact_exist_temp);
    lick_contact_exist_temp(empty_index) = {logical(0)};
    output = output & cell2mat(lick_contact_exist_temp);
end
end

function output = filter_sp_contact2_absent(sp_contact2_absent, lick_index, lick_nums)

% find trials where licks made contact with the spout
output = ones(1, numel(lick_index));
for j = 1:numel(lick_nums)
    sp_contact2_absent_exist_temp = cellfun(@(x, y) (y(x == lick_nums(j))), lick_index, sp_contact2_absent, 'UniformOutput', false);
    empty_index = cellfun('isempty', sp_contact2_absent_exist_temp);
    sp_contact2_absent_exist_temp(empty_index) = {logical(0)};
    output = output & cell2mat(sp_contact2_absent_exist_temp);
end
end

function output = filter_contact_centroid(contact_centroid, lick_index, lick_nums)

% find trials where contact centroid exists
output = ones(1, numel(lick_index));
for k = 1:numel(lick_nums)
    output = output & cell2mat(cellfun(@(x, y) numel([y{x == lick_nums(k)}]) > 1 & sum(isnan([y{x == lick_nums(k)}]), 'all') == 0, lick_index, contact_centroid, 'UniformOutput', false));
end
end

function output = filter_lick_exist(lick_index, lick_nums)

output = ones(1, numel(lick_index));
for k = 1:numel(lick_nums)
    output = output & cell2mat(cellfun(@(x) sum(x == lick_nums(k)) > 0, lick_index, 'UniformOutput', false));
end
end
