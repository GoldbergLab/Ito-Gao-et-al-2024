function videoData = loadVideoData(videoFilename, makeGrayscale, frames)
arguments
    videoFilename (1, :) char
    makeGrayscale (1, 1) logical = true
    frames double = []
end
%This currently works with grayscale avi and tif files

[~, ~, ext] = fileparts(videoFilename);

verbose = true;

try
    try
        % Check if user set NO_FFMPEG variable
        no_ffmpeg = evalin('base', 'NO_FFMPEG');
    catch ME
        % User did not set NO_FFMPEG environment variable
        no_ffmpeg = false;
    end
    if no_ffmpeg
        % User requested we not use ffmpeg, so skip this method
        if verbose
            disp('User requests no FFMPEG loading.')
        end
        throw(MException('User requests no FFMPEG'));
    end
    if verbose
        disp('Loading using fastVideoReader...');
    end
    videoData = fastVideoReader(videoFilename, [], frames);
catch
    if strcmp(ext, '.tif')
        % Check if file is a .tif file
        if verbose
            disp('Loading as .tif file')
        end
        tiffInfo = imfinfo(videoFilename);
        numFrames = length(tiffInfo);
        if isempty(frames)
            frames = 1:numFrames;
        end
        width = tiffInfo(1).Width;
        height = tiffInfo(1).Height;
        videoData = zeros([height, width, numFrames]);
        for k = frames
            videoData(:, :, k) = imread(videoFilename, k);
        end
    else
        try
            if verbose
                disp('Loading using read method with VideoReader')
            end
            video = VideoReader(videoFilename);
            if isempty(frames)
                videoData = read(video);
            else
                videoData = read(video, [min(frames), max(frames)]);
                videoData = videoData(:, :, frames-min(frames)+1);
            end
        catch
            try
                if verbose
                    disp('Loading using read method with VideoReader and native option')
                end
                video = VideoReader(videoFilename);
                if isempty(frames)
                    videoDataStruct = read(video, [1, video.NumFrames], 'native');
                else
                    videoDataStruct = read(video, [min(frames), max(frames)], 'native');
                end
                videoData = zeros([size(videoDataStruct(1).cdata), length(videoDataStruct)]);
                for k = 1:length(videoDataStruct)
                    videoData(:, :, k) = videoDataStruct(k).cdata;
                end
                if ~isempty(frames)
                    videoData = videoData(:, :, frames-min(frames)+1);
                end
            catch
                try
                    if verbose
                        disp('Loading using read method with VideoReader and native option with Inf as end frame')
                    end
                    video = VideoReader(videoFilename);
                    %    disp('Attempting to load video with method #1');
                    if isempty(frames)
                        videoDataStruct = read(video, [1, Inf], 'native');
                    else
                        videoDataStruct = read(video, [min(frames), Inf], 'native');
                    end
                    videoData = zeros([size(videoDataStruct(1).cdata), length(videoDataStruct)]);
                    for k = 1:length(videoDataStruct)
                        videoData(:, :, k) = videoDataStruct(k).cdata;
                    end
                    if ~isempty(frames)
                        videoData = videoData(:, :, frames-min(frames)+1);
                    end
                catch
                    if verbose
                        disp('Loading using VideoReader.readFrame and native option')
                    end
                    video = VideoReader(videoFilename);
                    videoData = zeros(video.Height, video.Width, int32(video.Duration * video.FrameRate));
                    frameNum = 1;
                    while hasFrame(video)
                        videoData(:, :, frameNum) = video.readFrame('native').cdata;
                        frameNum = frameNum + 1;
                        if ~isempty(frames) && frameNum > max(frames)
                            break
                        end
                    end
                    if ~isempty(frames)
                        videoData = videoData(:, :, frames);
                    end
                end
            end
        end
    end
end

% Get rid of duplicate RGB channels
videoDataSize = size(videoData);
if length(videoDataSize) == 4 && videoDataSize(3) == 1
    % For some reason we have a singleton dimension - let's get rid of it
    videoData = squeeze(videoData);
end
if length(videoDataSize) == 4 && videoDataSize(3) == 3 && makeGrayscale
    % third dimension is probably color channels
    videoData = squeeze(videoData(:, :, 1, :));
end