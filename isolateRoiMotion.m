function isolateRoiMotion(vid_path,v,data_path,Rotation,Centroid,imInvert)
% Create movies of the moving objects in a roi


%% Parameters

% Max number of frames to analyze for mean image
maxFrames = 100;

% Angular coordinates for circular roi
phi = linspace(0,2*pi,500);

% Rerun all steps in the analysis
force_redo = 0;


makeVid = 0;

%% Make mean image of roi & main image

if force_redo ||  isempty(dir([data_path filesep 'meanRoiImageData.mat']))
    
    % frames to use
    if length(Centroid.frames) > maxFrames
        framenums = Centroid.frames(1):...
            floor(length(Centroid.frames)/maxFrames): ...
            Centroid.frames(end);
    else
        framenums = Centroid.frames;
    end
    
    % Do not downsample the roi image
    dSample = 0;
    
    % Loop thru frames
    for i = 1:length(framenums)
        
        % Current frame
        cFrame = framenums(i);
        
        % Extract current parameters
        x = Centroid.x_pix(i);
        y = Centroid.y_pix(i);
        r = Centroid.r_pix;
        theta = Centroid.theta;
        tform = Rotation(i).tform_roi;
       
         % Current image
        im = getFrame(vid_path,v,cFrame,imInvert,'gray');
        
        % Extract roi
        [bw_mask,im_roi,roi_rect,bw_roi,imStable] = giveROI('circular',...
            im,x,y,r,theta,dSample,tform);
        
        % Add pixel values
        if i==1
            % Sum matrix for whole image
            imSum = double(im);
            
            % Sum matrix for stablized roi
            imSum_roi = double(imStable);
        else
            % Sum matrix for whole image
            imSum  = imSum + double(im);
            
            % Sum matrix for stablized roi
            imSum_roi  = imSum_roi + double(imStable);
        end
               
        % Clear for next
        clear x y r tform cFrame im theta im bw_mask im_roi roi_rect 
        clear bw_roi imStable
        
        % Update status
        disp(['Mean roi image: done ' num2str(i) ' of ' ...
            num2str(length(framenums))]);
    end
    
    % Calculate mean from sum image
    imMean = uint8(round(imSum./length(framenums)));
    
    % Calculate mean from sum image (stable roi)
    imRoiMean = uint8(round(imSum_roi./length(framenums)));
    
    % Save image image data (roi)
    save([data_path filesep 'meanRoiImageData'],'imRoiMean');
    
    % Save image image data
    save([data_path filesep 'meanImageData'],'imMean');
    
else
    % Load imRoiMean
    load([data_path filesep 'meanRoiImageData.mat'])
    
    % Load imMean
    load([data_path filesep 'meanImageData.mat'])
end


%% Mean image subtraction (roi)

if force_redo || isempty(dir([data_path filesep 'roiSequence.mat']))
    
    framenums = Centroid.frames;
    
    % Do not downsample the images
    dSample = 0;
    
    % Loop thru frames
    for i = 1:length(framenums)
        
        % Current frame
        cFrame = framenums(i);
        
        % Current parameters
        x = Centroid.x_pix(i);
        y = Centroid.y_pix(i);
        r = Centroid.r_pix;
        tform = Rotation(i).tform_roi;
        theta = Centroid.theta;
        
        % Current image
        im = getFrame(vid_path,v,cFrame,imInvert,'gray');
        
        % Extract roi
        [bw_mask,im_roi,roi_rect,bw_roi,imStable] = giveROI('circular',...
            im,x,y,r,theta,dSample,tform);

        % Adjust contrast 
        if i==1
            imRoiMean = imadjust(adapthisteq(imRoiMean));
        end
        imStable = imadjust(adapthisteq(imStable));
        
        % Subtract background
        warning off
        if imInvert
            im_seq(:,:,i) = imadjust(uint8(imcomplement(imsubtract(imRoiMean,imStable))));
        else
            im_seq(:,:,i) = imadjust(uint8(imcomplement(imsubtract(imStable,imRoiMean))));
        end
        warning on
        
        % Visualize
        if 0
            subplot(2,2,1)
            imshow(imStable,'InitialMagnification','fit');
            title('imStable')
            subplot(2,2,2)
            imshow(imRoiMean,'InitialMagnification','fit');
            title('mean image')
            subplot(2,2,3)
            warning off
            imshowpair(imStable,imRoiMean)
            warning on
            title('stable over mean image')
            subplot(2,2,4)
            imshow(im_seq(:,:,i),'InitialMagnification','fit');
            title('mean image subtraction')
            %title([num2str(cFrame)]
            pause(0.001);
        end
        
        % Clear for next
        clear x y r tform cFrame im theta im bw_mask im_roi roi_rect 
        clear bw_roi imStable
        
        % Update status
        disp(['Mean image subtraction: done ' num2str(i) ' of ' ...
            num2str(length(framenums))]);       
    end
    
    % Save image image data
    save([data_path filesep 'roiSequence'],'im_seq');
    
else
    % Load imRoiMean
    load([data_path filesep 'roiSequence.mat'])
end


%% Interactive mode

if force_redo || isempty(dir([data_path filesep 'blobParam.mat']))
    
    justReplay = 0;
    
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
            
            % Save params
            save([data_path filesep 'blobParam.mat'],'blobParam')
            
            % Break loop
            break
            
        elseif strcmp(choice,'Replay')
            justReplay = 1;
        end
    end
    
else
    % Load blobParam, parameter values
    load([data_path filesep 'blobParam.mat'])
end


%% Transform blobs back to global FOR

if isempty(dir([data_path filesep 'Blob data.mat']))
    
    framenums = Centroid.frames;
    
    for i = 1:length(Rotation)
        
        % Get roi image
        cIm = im_seq(:,:,i);
        
        % Find blobs in roi
        [props,bw_roi,areas,xB,yB] = findBlobs(cIm,blobParam.tVal,...
            'area',blobParam.areaMin,blobParam.areaMax);
        
        rMin = inf;
        rMax = 0;
        
        for j = 1:length(props)
            %rMin = min([rMin props(j).MinorAxisLength/2]);
            %rMax = max([rMax props(j).MajorAxisLength/2]);
            rVals(j,1) = props(j).MajorAxisLength/2;
        end
        
        rMax = ceil(quantile(rVals,0.9));
        rMin = floor(rMax/3);
        
        
        % Extract current parameters for whole image
        cFrame = framenums(i);
        x = Centroid.x_pix(i);
        y = Centroid.y_pix(i);
        r = Centroid.r_pix;
        theta = Centroid.theta;
        tform = Rotation(i).tform_roi;
        
        % Current whole frame
        im = getFrame(vid_path,v,cFrame,imInvert,'gray');
        
        % Get roi data
        [bw_mask,im_roi,roi_rect,bw_roi_mask] = giveROI('circular',im,x,y,r,theta,0);
        
        % Blobs in the G FOR
        bw_blobs_G = transCoord2d('bw L2G',tform,bw_roi,bw_mask,bw_roi_mask);
        
%         [cntrs,radii] = imfindcircles(bw_roi,[max([1 floor(rMin)]) ceil(rMax)],...
%                                       'ObjectPolarity','bright');
%           viscircles(cntrs, radii,'Color','b');  

        % Survey blobs
        propsG = regionprops(bw_blobs_G,'Centroid','Area',...
            'MajorAxisLength','MinorAxisLength',...
            'PixelIdxList','PixelList');
        
        % Store blob data
        B(i).propsG = propsG;
        B(i).propsL = props;
        
        if 0
            % Current whole frame
            imI = getFrame(vid_path,v,cFrame,~imInvert,'gray');
            
            h = imshow(imI,'InitialMag','fit');
            hold on
            
            % Make a truecolor all-green image, make non-blobs invisible
            green = cat(3, zeros(size(im)), ones(size(im)), zeros(size(im)));
            h = imshow(green,'InitialMag','fit');
            set(h, 'AlphaData', bw_blobs_G)
            
            title(['Frame ' num2str(framenums(i))]);
            
            pause(0.001)
            
        else
            disp(['Overlaying blobs: ' num2str(i) ' of ' num2str(length(Rotation))]);
        end
        
        clear props propsG bw_mask im_roi roi_rect bw_roi_mask bw_roi 
        clear areas xB yB
    end
    
    % Save image image data
    save([data_path filesep 'Blob data'],'B');
    
else
    % Load 'bwG_seq'
    load([data_path filesep 'Blob data']);
end

%% Filter out moving blobs

winLen = 10;

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

