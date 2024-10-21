function [cropped_mask, xlimits, ylimits, zlimits] = crop3Mask(mask, pad)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% crop3Mask: Crop the 3D mask to the bounding box containing true values
% usage:  cropped_mask = crop3Mask(mask)
%         [cropped_mask, xlimits, ylimits, zlimits] = crop3Mask(mask)
%
% where,
%    mask is a 3D logical array
%    cropped_mask is a 3D logical array containing only the smallest region 
%       of "mask" hat contains true values
%    xlimits is a 1x2 vector containing the minimum and maximum x
%       coordinates of true values within the mask.
%    ylimits is a 1x2 vector containing the minimum and maximum y
%       coordinates of true values within the mask.
%    zlimits is a 1x2 vector containing the minimum and maximum z
%       coordinates of true values within the mask.
%    pad is an optional integer indicating how much extra mask to take
%       around the edge of the cropped region. If the padding overlaps the 
%       edge of the mask, then less padding will be produced. Default is no
%       padding.
%
% <long description>
%
% See also: <related functions>
%
% Version: <version>
% Author:  Brian Kardon
% Email:   bmk27=cornell*org, brian*kardon=google*com
% Real_email = regexprep(Email,{'=','*'},{'@','.'})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('pad', 'var') || isempty(pad)
    pad = 0;
end

% Find the minimum and maximum coordinates for the true values along each
% axis.
[xlimits, ylimits, zlimits] = get3MaskLim(mask);

if isempty(xlimits) || isempty(ylimits) || isempty(zlimits)
    cropped_mask = [];
    return 
end

if pad ~= 0
    % Pad limits
    xlimits = xlimits + [-pad, pad];
    ylimits = ylimits + [-pad, pad];
    zlimits = zlimits + [-pad, pad];

    % Ensure limits do not exceed boundaries of mask
    xlimits(1) = max(1, xlimits(1));
    ylimits(1) = max(1, ylimits(1));
    zlimits(1) = max(1, zlimits(1));
    xlimits(2) = min(size(mask, 1), xlimits(2));
    ylimits(2) = min(size(mask, 2), ylimits(2));
    zlimits(2) = min(size(mask, 3), zlimits(2));
end

cropped_mask = mask(xlimits(1):xlimits(2), ylimits(1):ylimits(2), zlimits(1):zlimits(2));