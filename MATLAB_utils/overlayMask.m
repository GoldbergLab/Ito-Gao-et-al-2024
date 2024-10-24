function overlayImage = overlayMask(image, mask, color, transparency, origin)
% overlayMask: create an overlay image using an input image and a mask
% usage:  overlayImage = overlayMask(im, mask, color, transparency)
%
% where,
%    image is one of the following:
%       - A 2D H x W grayscale image array
%       - A 3D H x W x 3 color image array
%       - a 3D N x H x W stack of grayscale images
%       - a 4D N x H x W x 3 stack of RGB images
%       If it is a double array, pixel values should be in the range [0, 1]
%    mask is one of the following:
%       - A 2D binary mask of dimensions H x W
%       - A 3D binary mask stack of dimensions N x H x W
%       - A cell array of one of the above
%       If a mask stack is passed, image must also be a stack of the same
%       depth.
%    color is a color value to use to overlay the overlayMask, or a cell
%       array of colors, one per mask. Each color should be a 1x3 array of
%       RGB color values.
%    transparency is a transparency value from 0 to 1 to use to overlay 
%       the overlayMask, or a cell array 
%    origin is an optional 2-vector indicating where in the image the 
%       upper left pixel of the mask should be placed, or a cell array
%       of origin 2-vectors, one per mask. Default is [1, 1].
%
% overlayMask takes images and masks and overlays the masks using the
%   specified color and transparency. It can handle a variety of
%   combinations of 2D, 3D, and 4D images and masks.
%
% See also: <related functions>

% Version: <version>
% Author:  Brian Kardon
% Email:   bmk27=cornell*org, brian*kardon=google*com
% Real_email = regexprep(Email,{'=','*'},{'@','.'})

if ~exist('origin', 'var') || isempty(origin)
    origin = [1, 1];
end

if ~iscell(mask)
    mask = {mask};
end
if ~iscell(color)
    color = repmat({color}, 1, length(mask));
end
if ~iscell(origin)
    origin = repmat({origin}, 1, length(mask));
end
if ~iscell(transparency)
    transparency = repmat({transparency}, 1, length(mask));
end

% Check for mismatched #s of transparency, origin, and color
if length(transparency) ~= length(mask)
    error('If a cell array of transparency values is passed, the # of transparency values must match the # of masks.');
end
if length(origin) ~= length(mask)
    error('If a cell array of origins is passed, the # of origins must match the # of masks.');
end
if length(color) ~= length(mask)
    error('If a cell array of colors is passed, the # of colors must match the # of masks.');
end

if any(size(image) == 1) 
    image = squeeze(image);
end

for k = 1:length(mask)
    if any(size(mask{k}) == 1)
        mask{k} = squeeze(mask{k});
    end
end

for k = 1:length(color)
    % Check that colors are valid.
    if length(color{k}) ~= 3 || ~isnumeric(color{k})
        error('Expected colors to each be a 1x3 array of RGB color values.')
    end
end

switch ndims(image)
    case 2
        nImageFrames = 1;
        nImageChannels = 1;
    case 3
        if size(image, 3) == 3
            % Probably a single RGB image
            nImageFrames = 1;
            nImageChannels = 3;
        else
            % Must be a stack of grayscale images
            nImageFrames = size(image, 1);
            nImageChannels = 1;
        end
    case 4
        if size(image, 4) == 3
            nImageFrames = size(image, 1);
            nImageChannels = 3;
        else
            error('Invalid image size: %s', num2str(size(image)));
        end
    otherwise
        error('Invalid image size: %s', num2str(size(image)));
end
if nImageChannels == 1
    % Convert grayscale image to color
    image = cat(length(size(image))+1, image, image, image);
end

% Check that transparency is in [0, 1]
for k = 1:length(transparency)
    if transparency{k} < 0 || transparency{k} > 1
        error(['transparency must be in the range [0, 1]. Instead, transparency was ', transparency{k}])
    end
end

for k = 1:length(mask)
    if isempty(mask{k})
        continue;
    end
    if ~islogical(mask{k})
        mask{k} = logical(mask{k});
    end
%     [mask{k}, ylimits, xlimits] = cropMask(mask{k});
%     if isempty(mask{k})
%         continue;
%     end
%     if isrow(color{k})
%         color{k} = color{k}';
%     end
%     origin{k} = origin{k} + [xlimits(1), ylimits(1)] - 1;
    switch ndims(mask{k})
        case 2
            if nImageFrames == 1
                % Image is a 3D color image
                % Mask is a single 2D image. Make it a 3D color image to
                % match image.
                mask{k} = tensorprod(double(mask{k}), color{k} * (1 - transparency{k}));
            else
                % Image is a 4D color stack
                % Mask is a single 2D image. Make it a 4D color stack to
                % match image.
                mask{k} = repmat(reshape(tensorprod(double(mask{k}), color{k} * (1 - transparency{k})), [1, size(mask{k}), 3]), [nImageFrames, 1, 1, 1]);
%                 mask{k} = permute(repmat(mask{k}, [1, 1, 3, nImageFrames]), [4, 1, 2, 3]);
            end
        case 3
            % Mask is a 3D stack. Make it a 4D color image
%             mask{k} = cat(4, mask{k}*color{k}(1), mask{k}*color{k}(2), mask{k}*color{k}(3));
            mask{k} = tensorprod(double(mask{k}), color{k} * (1 - transparency{k}));
%             mask{k} = repmat(mask{k}, [1, 1, 1, 3]);
        otherwise
            error('Invalid mask dimensions: %s', num2str(size(mask{k})))
    end
end

% Convert color input to an RGB triplet
%color = rgb(color);

% Transform overlayMask to an appropriate data type with appropriate
% scaling for the image

switch class(image)
    case {'single', 'double'}
        maxVal = max(image(:));
        minVal = min(image(:));
        if minVal < 0 || maxVal > 1
            error(['Input image is of type double - its values must be in the range [0, 1]. Instead it had the range [', num2str(minVal), ', ', num2str(maxVal), ']'])
        end
    case {'int8', 'int16', 'int32', 'int64', 'uint8', 'uint16', 'uint32', 'uint64'}
        image = double(image)/double(intmax(class(image)));
    otherwise
        error(['image must be a numeric class. Instead, it was', class(image)]);
end

overlayImage = image;
for k = 1:length(mask)
    if isempty(mask{k})
        continue;
    end
    thisOrigin = origin{k};

    w = size(image, length(size(image))-1);
    h = size(image, length(size(image))-2);
    wm = size(mask{k}, length(size(mask{k}))-1);
    hm = size(mask{k}, length(size(mask{k}))-2);

    [x1, x2, y1, y2, overlaps] = getRectangleOverlap(1, w, 1, h, thisOrigin(1), thisOrigin(1) + wm - 1, thisOrigin(2), thisOrigin(2) + hm - 1);
    if overlaps
        if (x1 == 1 && x2 == w && y1 == 1 && y2 == h)
            overlayImage = overlayImage + mask{k};
            disp('simple overlay')
        else
            x1m = x1 - thisOrigin(1) + 1;
            y1m = y1 - thisOrigin(2) + 1;
            x2m = x2 - thisOrigin(1) + 1;
            y2m = y2 - thisOrigin(2) + 1;
            if nImageFrames == 1
                overlayImage(y1:y2, x1:x2, :) = overlayImage(y1:y2, x1:x2, :) + mask{k}(y1m:y2m, x1m:x2m, :);
            else
                overlayImage(:, y1:y2, x1:x2, :) = overlayImage(:, y1:y2, x1:x2, :) + mask{k}(:, y1m:y2m, x1m:x2m, :);
            end
        end
    end
end
