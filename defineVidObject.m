function v = defineVidObject(vid_path,ext)
% Identifies the properties of the video and returns them in a structure
%    vid_path - path to the video file or directory of image files
%    ext - file extension (e.g. 'jpg')
%
% Note: currently only handles jpgs for image sequences
%
% Developed by McHenryLab at UC Irvine

% Extension default
if nargin < 2
    ext = 'jpg';
end




% If directory of images . . .
if isdir(vid_path)
     
    % Get image sequence
    a = dir([vid_path filesep '*.' ext]);
    
    % Check for files
    if isempty(a)
        error(['No images found with the ' ext ' extension in ' ...
              vid_path]);
    end
    

    % Read one image for dimensions
    im = imread([vid_path filesep a(1).name]);
    
    % Fill in basics for video info
    v.Path   = vid_path;
    v.Width  = size(im,2);
    v.Height = size(im,1);
    
    
    % Loop thru files
    for i = 1:length(a)
        
        % Get filename
        fName = a(i).name;
        
        % Get frame number
        frNum(i,1) = str2num(fName((find(fName=='_')+1):(end-length(ext)-1)));
        
    end
    
    % Store away number of frames and info on files
    v.UserData.FirstFrame = min(frNum);
    v.UserData.LastFrame = max(frNum);
    v.UserData.FileInfo  = a;
    
% If single file . . .
else
    
    % Extract file parts
    %[pathstr,name,ext] = fileparts(vid_path);
    
    % Find video info
    v = VideoReader(vid_path);
    
    % Store number of frames
    v.UserData.FirstFrame = 1;
    v.UserData.LastFrame = floor(v.FrameRate*v.Duration);
end



