function imMean = makeMeanImage(vid_path,newMean)
% makeMeanImage generates a mean background image for background
% subtraction. Two options are available: mean image based on displacement
% of objects in video or static background with no objects
%
% INPUTS: 

%% Parameters

% Max number of frames for creating the mean image
maxFrames = 500;

% If directory of images . . .
if isdir(vid_path)
    
    % Get image sequence
    a = dir([vid_path filesep '*.jpg']);
    
    % total number of frames in 'a'
    frTot = length(a);
else
    % Give error message
    error('Video path not given')
end
    
%% Create or load mean image
    
% Look for mean image
a2 = dir([vid_path filesep 'meanImage.tif']);

% Calculate mean image if it does not exist & we want to use meanImage
if isempty(a2) && ~(newMean)
    
    % Define list of frame numbers, depending on max number of frames
    % requested
    if frTot > maxFrames
        dframe = floor(frTot/maxFrames);
        frIdx = 1:dframe:frTot;
        clear dframe
    else
        frIdx = 1:frTot;
    end
    
    % Create waitbar
    h = waitbar(0,['Mean image: ' num2str(1)],'Name','Mean image',...
                     'CreateCancelBtn',...
                     'setappdata(gcbf,''canceling'',1)');
        setappdata(h,'canceling',0)
    
    % Create sum image based on first frame
    [imCurr,tmp] = imread([vid_path filesep a(1).name]);
    
    % Start sum image
    imSum = double(imCurr.*0);
    
    % Min image (hilights?) 
    if strcmp(class(imCurr),'uint8')
        imMin = imCurr.*0 + 255; 
        
    elseif strcmp(class(imCurr),'uint16')
        imMin = imCurr.*0 + 2.^16;
        
    else
        error('Code only supports 16 and 8 bit images')
        
    end
    
    clear tmp
      
    % Loop through frames (iteratively updates 'imSum' and 'imMin') 
    for i = 1:length(frIdx)
        
        % Add current frame to sum image
        [imCurr,tmp] = imread([vid_path  filesep a(i).name]);
        
        % Sum pixel values (cumulative sum)
        imSum  = imSum + double(imCurr);
       
        % Update min image
        tMin(:,:,1) = imMin;
        tMin(:,:,2) = imCurr;
        imMin = min(tMin,[],3);
        
        % Clear for next
        clear tmp 
        
        % Check for Cancel button press
        if getappdata(h,'canceling')
            close(h,'force')
            error('Execution stopped by user');
            
            % Otherwise, update status
        else
            waitbar(i/length(frIdx),h,...
                ['Mean image: ' num2str(i) ' of ' ...
                num2str(length(frIdx)) ' frames'])
        end       
    end
    
    if strcmp(class(imCurr),'uint8')
        
        % Calculate mean from sum image
        imMean = uint8(round(imSum./length(frIdx)));
    elseif strcmp(class(imCurr),'uint16')
        
        % Calculate mean from sum image
        imMean = uint16(round(imSum./length(frIdx)));
    end
    
    % Write mean image to movie dir
    imwrite(imMean,[vid_path filesep 'meanImage.tif'],'tif','Compression','none');

    % Close status bar
    close(h,'force')

% create a new background image, and use as mean image
elseif newMean
    
    % Look for meanImage2
    a3 = dir([vid_path filesep 'meanImage2.tif']);
    
    % Create meanImage2 if it does not already exist
    if isempty(a3)
        
        % Read first frame of image sequence
        im = imread([vid_path filesep a(1).name]);
        
        % Select region of interest around pred
        warning off
        imshow(im)
        warning on
        title('Choose ROI around moving objects (e.g. fish)')
        
        % Interactively find ROI
        h = impoly;
        roi_poly = wait(h);
        
        % Store results
        tmp = getPosition(h);
        roi.x = tmp(:,1);
        roi.y = tmp(:,2);
        
        delete(h), close all;
        
        % create a binary mask based on the selected ROI (select fish)
        maskFish = roipoly(im,roi.x,roi.y);
        
        % estimate background image
        imMean = roifill(im,maskFish);
        
        % Write mean image to movie dir
        imwrite(imMean,[vid_path filesep 'meanImage2.tif'],'tif','Compression','none');
        
        % ...otherwise load meanImage2
    else
        imMean = imread([vid_path filesep 'meanImage2.tif']);
    end
    
    disp('   Using meanImage2 for bkgnd subtraction...')
    

% Load images, if present
else

    imMean = imread([vid_path filesep 'meanImage.tif']);  

end




