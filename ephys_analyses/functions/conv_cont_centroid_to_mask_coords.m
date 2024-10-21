function contact_centroid_mask = conv_cont_centroid_to_mask_coords(contact_centroid, fiducial_xy)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function converts contact_centroid into mask coordinates.
%
% contact_centroid is an array containing contact_centroid, extracted from
% sp_struct, where the first dimension should be neuron and the second
% dimension is trials. 

% Author:  Brendan Ito
% Email:   bsi8@cornell.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define fiducial and transform fiducial_x into mask coordinates
fiducial_xy(:, 1) = fiducial_xy(:, 1) - (400 - 240);

% convert contact centroid coordinates to mask coordinates
for i = 1:size(contact_centroid, 1)
    for j = 1:size(contact_centroid, 2)
        contact_centroid_mask{i, j}(:, 1) = contact_centroid{i, j}(:, 1) - fiducial_xy(1);
        contact_centroid_mask{i, j}(:, 2) = contact_centroid{i, j}(:, 2) - fiducial_xy(2);
        contact_centroid_mask{i, j}(:, 3) = contact_centroid{i, j}(:, 3);
    end
end