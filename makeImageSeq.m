function makeImageSeq(filename,pathname,varargin)
% Creates a grayscale image sequence from a video file
%   filename - name of video file
%   pathname - path to video file

% makeImageSeq(filename,pathname) 
% writes image sequnence to pathname

% makeImageSeq(filename,pathname,newpath) 
% writes image sequnence to newpath

%
% Developed by McHenryLab at UC Irvine

% Prompt to browse to filename
if nargin < 1
    [filename,pathname] = uigetfile({'*.MP4';'*.mp4';'*.MOV'});
end

if nargin>2
    newpath = varargin{1};
else
    newpath = pathname;
end

% Get file info
[pathstr,name,ext] = fileparts([pathname filesep filename]);

% Get video information
v = defineVidObject([pathname filesep filename]);

% Make directory with same name
mkdir(pathname,name);

% Save video info
save([pathname filesep name filesep 'vid_info'],'v');

% Define grayscale colormap
cmap = colormap('gray');

% loop thru frames
for i = 1:floor(v.FrameRate * v.Duration)

    % Read next available frame
    im = readFrame(v,'native');
    
    % Convert to grayscale, enhance contrast
    imGray = rgb2gray(im);
    
    % Frame number string
    frnum = ['0000000' num2str(i)];
    frnum = frnum(end-6:end);
    
    % Create filename for frame
    newname = ['frame_' frnum '.jpg'];
    
    % Write image
    %imwrite(imGray,cmap,[newpath filesep name filesep newname],'jpg','BitDepth',8);
    imwrite(imGray,[newpath filesep name filesep newname],'jpg','BitDepth',8);
    
    % Update status
    disp(['  Written ' num2str(i) ' of ' ...
          num2str(floor(v.FrameRate * v.Duration)) ' frames'])
end