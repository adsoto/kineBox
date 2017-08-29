function v = defineVidObject(vid_path)
% Identifies the properties of the video and returns them in a structure
%    vid_path - path to the video file or directory of image files
%
% Note: currently only handles jpgs for image sequences
%
% Developed by McHenryLab at UC Irvine

% If directory of images . . .
if isdir(vid_path)
    
    % Get image sequence
    a = dir([vid_path filesep '*.jpg']);
    
    % If there is a video info dir . . .
    if ~isempty(dir([vid_path filesep 'vid_info.mat']))
        % Load video info (v)
        load([vid_path filesep 'vid_info'])
    else
        % Read one image for dimensions
        im = imread([vid_path filesep a(1).name]);
        
        % Fill in basics for video info
        v.Path   = vid_path;
        v.Width  = size(im,2);
        v.Height = size(im,1); 
    end    
    
    % Store away number of frames and info on files
    v.UserData.NumFrames = length(a);
    v.UserData.FileInfo  = a;
    
% If single file . . .
else
    % Find video info
    v = VideoReader(vid_path);
    
    % Store number of frames
    v.UserData.NumFrames = floor(v.FrameRate*v.Duration);
end



