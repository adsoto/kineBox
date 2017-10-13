function im = getFrame(vid_path,v,fr_num,imInvert,clrMode,imMean)  
% Reads image file from a video. 
%   v - video structure
%   fr_num - frame number desired
%   vid_path - path to video (use path in v structure as default)
%   imInvert - logical that indicates whether to invert a grayscale image
%   clrMode - 'rgb' or 'gray' format for image 
%   imMean - (optional) mean image of video sequence

%
% Developed by McHenryLab at UC Irvine

% If it is an image sequence . . .
if isfield(v.UserData,'FileInfo');
    
    % Default for image invert
    if nargin < 4
        imInvert = 0;
    end
    
    % Default color mode
    if nargin < 5
        clrMode = 'rgb';
    end
    
    if nargin < 6
        imMean = [];
    end
    
    % If no frame number given . . .
    if nargin < 3 || isempty(fr_num)
        warning('No frame number given, defaulting to first frame')
        fr_num = 1;
        
    % Otherwise . . .
    else
        % Check requested frame number
        if fr_num>v.UserData.LastFrame
            error(['Video sequence does not have a frame ' num2str(fr_num)]);
        end
        
        % Frame index
        iFrame = fr_num-v.UserData.FirstFrame + 1;
        %iFrame = fr_num;
        
        % Get filename and extension
        fName = v.UserData.FileInfo(iFrame).name;
        ext   = fName(find(fName=='.')+1:end);
        fName = fName(1:(find(fName=='.')-1));
       
        % Get frame number
        frNum = str2num(fName((find(fName=='_')+1):end));
        
        % Check for match
        if frNum~=fr_num
            error('file numbering does not match the video info');
        end
    end
    
    % Read image
    im = imread([vid_path filesep v.UserData.FileInfo(iFrame).name]);    
    
    % Convert to grayscale
    if strcmp(clrMode,'gray')
        im = rgb2gray(im);
    end
    
    % Subtract mean image, if present
    if ~isempty(imMean)
        im = imsubtract(im,imMean);
    end
    
    % Invert, if requested
    if imInvert
        im = imcomplement(im);
    end
  
% If a video file . . .
else
    if nargin > 2
        warning('Cannot specify image number for a video.  Reading next image.')
    end
    % Read next available frame
    im = readFrame(v);
end
