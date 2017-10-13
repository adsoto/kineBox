function isolateRoiMotion(vid_path,v,data_path,Rotation,Centroid,imInvert)
% Create movies of the moving objects in a roi


%% Parameters

% Max number of frames to analyze for mean image
maxFrames = 100;

% Angular coordinates for circular roi
phi = linspace(0,2*pi,500);

% Rerun all steps in the analysis
force_redo = 0;

% Make movie of resulting data 
makeVid = 0;


%% Make mean images of roi & main image

if 1 %force_redo ||  ...
        (isempty(dir([data_path filesep 'meanRoiImageData.mat'])) && ...
         isempty(dir([data_path filesep 'meanImageData.mat'])) )
    
     % Make mean images 
    [imMean,imRoiMean] = motionImage(vid_path,v,'mean roi','none',...
                                  maxFrames,Centroid,Rotation,0);
  
    % Plot mean images
    f = figure;
    subplot(3,2,[1 3])
    imshow(imMean,'InitialMag','fit');
    title('Mean image')
    subplot(3,2,5)
    imhist(imMean)
    
    subplot(3,2,[2 4])
    imshow(imRoiMean,'InitialMag','fit');
    title('Mean roi image')
    subplot(3,2,6)
    imhist(imRoiMean)
    
    % Prompt to proceed
    ButtonName = questdlg('Do the mean images look good?', ...
        '', 'Yes, proceed', 'No', 'Cancel', 'Yes, proceed');
    
    switch ButtonName,
        case 'Cancel',
            disp(''); disp('Quitting analysis . . .');beep
            return
            
        case 'No',
            disp(''); disp('Hmm take a look at the code to figure out what went wrong.');
            beep
            return
    end % switch
    
    close(f)    
    clear ButtonName
    
    % Save image image data (roi)
    save([data_path filesep 'meanRoiImageData'],'imRoiMean');
    
    % Save image image data
    save([data_path filesep 'meanImageData'],'imMean');
    
else
    disp('Loading mean images . . .')
    
    % Load imRoiMean
    load([data_path filesep 'meanRoiImageData.mat'])
    
    % Load imMean
    load([data_path filesep 'meanImageData.mat'])
end

disp(' ');


%% Mean image subtraction (roi)

if 1 %force_redo || isempty(dir([data_path filesep 'roiSequence.mat']))
    
    % Get sequence of roi with subtracted mean images
    im_seq = giveSequence(vid_path,v,'mean roi sub','enhance contrast',...
        imRoiMean,Centroid.frames,Centroid,Rotation,0,1);

    % Prompt to proceed
    ButtonName = questdlg('Ready to play sequence?', ...
        '', 'Yes','Cancel', 'Yes');
    
    switch ButtonName,
        case 'Cancel',
            disp(''); disp('Quitting analysis . . .');beep
            return
    end % switch
    clear ButtonName
      
    f = figure;
    
    while true
        
        % Loop thru frames
        for i = 1:size(im_seq,3)
            
            % Display image
            imshow(im_seq(:,:,i),'InitialMag','fit');
            
            % Pause to render
            pause(0.1)           
        end
        
        % Prompt to proceed
        ButtonName = questdlg('Proceed to next step?', ...
            '', 'Yes','Replay','Cancel', 'Yes');
        
        switch ButtonName,
            case 'Yes'
                close(f)
                break
                
            case 'Cancel',
                disp(''); disp('Quitting analysis . . .');beep
                return
        end % switch
        clear ButtonName
        
    end
    
    % Save image image data
    save([data_path filesep 'roiSequence'],'im_seq');
    
else
    
    disp('Loading im_seq . . .')
    
    % Load im_seq
    load([data_path filesep 'roiSequence.mat'])
end

disp(' ');


%% Interactive mode: select threshold and area limits

if force_redo || isempty(dir([data_path filesep 'blobParam.mat']))
    
    justReplay = 0;
    
    f = figure;
    
    while true
        
        if ~justReplay
            i = 1;
            iMax = min([50 size(im_seq,3)]);
            
            % Find threshold
            tVal = imInteract(im_seq(:,:,i),'threshold');
            
            % Find area
            [areaMin,areaMax] = imInteract(im_seq(:,:,i),'area',tVal);
        end
        
        % Play
        for i = 1:iMax
            imInteract(im_seq(:,:,i),'display blobs',tVal,areaMin,areaMax);
            title(['Frame ' num2str(Centroid.frames(i))]);
            pause(0.1);
        end
        
        % Reset
        justReplay = 0;
        
        % Prompt
        choice = questdlg('Do threshold and area values look good?', ...
            '', 'Yes, proceed','No, redo','Replay','Yes, proceed');
        
        % Break, if good
        if strcmp(choice,'Yes, proceed')
            
            % Store parameters
            blobParam.tVal = tVal;
            blobParam.areaMin = areaMin;
            blobParam.areaMax = areaMax;
            
            % Close figure window
            close(f)
            
            % Save params
            save([data_path filesep 'blobParam.mat'],'blobParam')
            
            % Break loop
            break
            
        elseif strcmp(choice,'Replay')
            justReplay = 1;
        end
    end
    
else
    disp('Loading blobParam . . .');
    
    % Load blobParam, parameter values
    load([data_path filesep 'blobParam.mat'])
end

disp(' ');


%% Transform blobs back to global FOR

if force_redo || isempty(dir([data_path filesep 'Blob data.mat']))
    
    % Get blobs
    B = anaBlobs(vid_path,v,'G&L props',...
        im_seq,Centroid.frames,Centroid,Rotation,blobParam,imInvert);
    
    % Save image image data
    save([data_path filesep 'Blob data'],'B');
    
else
    
    disp('Loading B structure . . .');
    
    % Load 'B' structure
    load([data_path filesep 'Blob data']);
end

disp(' ');


%% Filter out moving blobs

%TODO: Fix the code below . . .

if 1 %force_redo || isempty(dir([data_path filesep 'Blob, filtered data.mat']))
    
    winLen = 10;
    
    %tVal = imInteract(imcomplement(imStackScore),'threshold');
    tVal = 0.5871;
    
    % Find blobs that don't move in global FOR
    Bf = anaBlobs(vid_path,v,'filter motion',B,winLen,Centroid,...
        Rotation,blobParam,tVal);
    
    aniData(vid_path,v,'blobs G&L',Bf)
    
    % Save blob data
    save([data_path filesep 'Blob, filtered data'],'Bf');
    
else
     disp('Loading Bf structure . . .');
    
    % Load 'Bf' structure
    load([data_path filesep 'Blob, filtered data']);
end

disp(' ');disp(' ')

return

bw_seq = giveSequence(vid_path,v,'filter blob motion',im_seq,B,fr_num,...
    winLen,'none',Centroid,Rotation,dSample,imInvert)

framenums = Centroid.frames;

for i = ceil(winLen/2):length(Rotation)-ceil(winLen/2)
    
    % Get roi image
    cIm = im_seq(:,:,i);

    % Extract current parameters for whole image
    cFrame = framenums(i);
    x = Centroid.x_pix(i);
    y = Centroid.y_pix(i);
    r = Centroid.r_pix;
    theta = Centroid.theta;
    tform = Rotation(i).tform_roi;
    
    % Current whole frame
    im = getFrame(vid_path,v,cFrame,imInvert,'gray');
    
    startFrame = max([1 i-floor(winLen/2)]);
    endFrame = min([length(Rotation) i+floor(winLen/2)]);
    
    winFrames = startFrame:endFrame;
    
    % Loop thru window of images
    for j = 1:length(winFrames)
        
        % Current frame
        iFrame = winFrames(j);
        
        % Start with blank
        currIm  = logical(zeros(size(im)));
        
        % Score pixels with blobs
        for k = 1:length(B(iFrame).propsG),
            currIm(B(iFrame).propsG(k).PixelIdxList) = 1;            
        end
        
        % Store resulting image
        bwStack(:,:,j) = currIm;
    end
    
    % image of the pixels that are most static
    imStackScore = uint8(sum(double(bwStack),3)./winLen.*255);
    
    %tVal = imInteract(imcomplement(imStackScore),'threshold');
    tVal = 0.5871;
    
    bwNew = im2bw(imStackScore,tVal);
    
    % Fill holes
    bwNew = imfill(bwNew,'holes');
    
    %bwPer = bwperim(bwNew);
    
%     % Get roi data
%     [bw_mask,im_roi,roi_rect,bw_roi_mask,imStable] = giveROI('circular',...
%         imStackScore,x,y,r,theta,0,tform);
    
    
    %imStackScore = uint8(mean(bwStack,3)./winLen.*255);
    
    %imshow(imStackScore)
    
    if 1
        % Current whole frame
         imI = getFrame(vid_path,v,cFrame,~imInvert,'gray');
        
        h = imshow(imI,'InitialMag','fit');
        hold on
%         
        % Make a truecolor all-green image, make non-blobs invisible
        green = cat(3, zeros(size(im)), ones(size(im)), zeros(size(im)));
        h = imshow(green,'InitialMag','fit');
        set(h, 'AlphaData', double(bwNew).*0.2)
        
        title(['Frame ' num2str(framenums(i))]);
        hold off
        
        % Log frame
        if makeVid
            M(i) = getframe(gcf);
        end
        
        
        pause(0.001)
        
    else
        disp(['Overlaying blobs: ' num2str(i) ' of ' num2str(length(Rotation))]);
    end
    
    clear props propsG bw_mask im_roi roi_rect bw_roi_mask bw_roi
    clear areas xB yB
end



 % Write movie to disk
    if makeVid
        vid_save_path = uigetdir(data_path,'Save movie');
        vInfo = VideoWriter([vid_save_path filesep 'video tubefeet.mp4'],'MPEG-4');
        vInfo.FrameRate = 15;
        open(vInfo)
        writeVideo(vInfo,M)
        close(vInfo)
    end

% function [im,roi_mask] = isolate_roi(im,x,y,r,theta)
% 
% % Maximum size of an image dimension
% maxSize = 250;
% 
% % rectangular ROI vector
% roi_rect = ceil([x-r y-r r*2 r*2]);
% 
% % Crop image
% im  = imcrop(im,roi_rect);
% 
% % Circular coordinates for new roi
% xC    = r.*cos(theta)+ size(im,2)/2;
% yC    = r.*sin(theta) + size(im,1)/2;
% 
% % White out pixels outside of circular roi
% roi_mask = roipoly(size(im,1),size(im,2),round(yC),round(xC));
% 
% % White out area around roi
% im(~roi_mask) = 255;
% 
% % Reduce size (helps speed up image registration)
% if length(im)>maxSize
%     
%     imFactor = maxSize/length(im);
%     
%     im = imresize(im,imFactor);
%     
%     % White out pixels outside of circular roi
%     roi_mask = roipoly(size(im,1),size(im,2),round(yC*imFactor),round(xC*imFactor));
% end
% 
% % Check for square
% if size(im,1)~=size(im,2)
%     error('Image not square');
% end


function meanRoiImage

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

function aniInteract(roi)

% Frame number
i = 1;

% First frame
im = roi(:,:,i);

% Plot image
figure
imshow(im,'InitialMagnification','fit');
title('Choose threshold w/arrows, select blob, press return')
hold on

% Guess threshold
tVal = graythresh(im);

% Loop interaction
while true
    
    % Make binary
    bw = im2bw(im,tVal);
    
    % Fill holes
    bw = imfill(bw,'holes');
    
    % Trace perimeter
    [y, x] = find(bwperim(bw,8));
    
    % Overlay
    h = plot(x,y,'.g');
    
    % Get input
    [x,y,b] = ginput(1);
    
    % Up arrow
    if b==30
        tVal = min([tVal+0.05 1]);
        
        % Down arrow
    elseif b==31
        tVal = max([tVal-0.05 0]);
        
  
    % Return
    elseif isempty(b)
        
        % Check for coordinates
        if isempty(xPos)
            warning('You need to select a blob before pressing return')
        else
            break
        end
    end
    
    % Remove perimeter
    delete(h)
end

props = regionprops(bw,'Centroid','Area',...
    'MajorAxisLength','MinorAxisLength');

