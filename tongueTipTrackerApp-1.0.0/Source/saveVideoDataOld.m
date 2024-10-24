function saveVideoData(videoData, filename, varargin)

if nargin == 3
    videoType = varargin{1};
else
    if ndims(videoData) == 3
        videoType = 'Grayscale AVI';
    elseif ndims(videoData) == 4
        videoType = 'Uncompressed AVI';
    end
end

v = VideoWriter(filename, videoType);
open(v);
videoSize = size(videoData);
% for k = 1:videoSize(3)
%     writeVideo(v, videoData(:, :, k));
% end
if isa(videoData, 'double') && find(videoData > 1, 1)
    videoData = uint8(videoData);
end
if ndims(videoData) == 3
    videoSizeAdjusted = [videoSize(1:2), 1, videoSize(3)];
    videoData = reshape(videoData, videoSizeAdjusted);
end
writeVideo(v, videoData);
close(v);