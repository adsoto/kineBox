function imMean = makeMeanImage(vid_path,v,preprocess,fr_nums)
% Creates a mean image from a video
%  v - structure of info about video
%
% Developed by McHenryLab at UC Irvine

% Whether to include waitbar for status
useWaitBar = 0;

% Check for image sequence
if ~isdir(vid_path)
    error('This function requires that the video is saved as a series of images')
end

% Set default for preprocessing
if nargin < 3
    preprocess = 'none';
end
   
% If frame numbers not provided . . .
if nargin < 4
    
    % Max number of frames to analyze
    maxFrames = 200;
    
    % Proportion of video duration to exclude, from the beginning
    exclude_prop = 0.2;
    
    % Full roi
    p.roi_x = [1 v.Width v.Width 1 1];
    p.roi_y = [1 v.Height 1 1 v.Height];
    
    % Define list of frame numbers, depending on max number of frames
    % requested
    if (1-exclude_prop)*v.UserData.NumFrames > maxFrames
        dframe = floor(v.UserData.NumFrames/maxFrames);
        frame1 = round(exclude_prop*v.UserData.NumFrames);
        fr_nums = frame1:dframe:v.UserData.NumFrames;
        clear dframe frame1
    else
        fr_nums = 1:v.UserData.NumFrames;
    end
end

if useWaitBar
    % Create waitbar
    h = waitbar(0,...
        ['Mean image: ' num2str(1)],...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
end

% Create sum image based on first frame
imCurr = getFrame(vid_path,v,1);  

% Convert to grayscale
imCurr = rgb2gray(imCurr);

% Enhance contrast, if requested
if strcmp(preprocess,'enhance contrast');
    imCurr = adapthisteq(imCurr);
end


imSum = double(imCurr);
clear imCurr 

% Loop through frames
for i = 1:length(fr_nums)
    
    cFrame = fr_nums(i);
    
    % Get current frame
    imCurr       = getFrame(vid_path,v,cFrame);  
    
    % Convert to grayscale
    imCurr = rgb2gray(imCurr);

    % Enhance contrast, if requested
    if strcmp(preprocess,'enhance contrast');
        imCurr = imadjust(imCurr);
    end

    imSum  = imSum + double(imCurr);
    clear imCurr
    
    if useWaitBar
        %Update status bar
        h = waitbar(i/length(frIdx),h,...
            ['Mean image: ' num2str(i) ' of ' num2str(length(frIdx)) ' frames']);
        
        % Quit m-file, if cancel button pushed
        if getappdata(h,'canceling')
            close force
            return
        end
    end
    
end

% Calculate mean from sum image
imMean = uint8(round(imSum./length(fr_nums)));

%imMean = imMean(:,:,1);

if useWaitBar
    close(h)
end

end



