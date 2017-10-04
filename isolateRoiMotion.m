function isolateRoiMotion(vid_path,v,data_path,Rotation,Centroid)
% Create movies of the moving objects in a roi


%% Parameters

% Max number of frames to analyze for mean image
maxFrames = 100;

% Angular coordinates for circular roi
phi = linspace(0,2*pi,500);


%% Make mean image of roi

if isempty(dir([data_path filesep 'meanRoiImageData.mat']))
    
    % frames to use
    if length(Centroid.frames) > maxFrames
        framenums = Centroid.frames(1):...
            floor(length(Centroid.frames)/maxFrames): ...
            Centroid.frames(end);
    else
        framenums = Centroid.frames;
    end
    
    % Loop thru frames
    for i = 1:length(framenums)
        
        % Current frame
        cFrame = framenums(i);
        
        % Extract current parameters
        x = Centroid.x_pix(i);
        y = Centroid.y_pix(i);
        r = Centroid.r_pix;
        tform = Rotation(i).tform_roi;
        
        % Load current image
        im = rgb2gray(getFrame(vid_path,v,cFrame));
        
        % Get roi image
        [im_roi,imMask] = isolate_roi(im,x,y,r,phi);
        
        % Stablize rotation
        im_roi = imwarp(im_roi,tform,'OutputView',imref2d(size(im_roi)));
        
        % White out beyond roi
        im_roi(~imMask) = 255;
    
        % Add pixel values
        if i==1
            imSum = double(im_roi);
        else
            imSum  = imSum + double(im_roi);
        end
        
        % Clear for next
        clear x y r tform cFrame im
        
        % Update status
        disp(['Mean roi image: done ' num2str(i) ' of ' ...
            num2str(length(framenums))]);
    end
    
    % Calculate mean from sum image
    imRoiMean = uint8(round(imSum./length(framenums)));
    
    % Save image image data
    save([data_path filesep 'meanRoiImageData'],'imRoiMean');
    
else
    % Load imRoiMean
    load([data_path filesep 'meanRoiImageData.mat'])
end


%% Mean image subtraction


if isempty(dir([data_path filesep 'roiSequence.mat']))
    
    framenums = Centroid.frames;
    
    % Loop thru frames
    for i = 1:length(framenums)
        
        % Current frame
        cFrame = framenums(i);
        
        % Current parameters
        x = Centroid.x_pix(i);
        y = Centroid.y_pix(i);
        r = Centroid.r_pix;
        tform = Rotation(i).tform_roi;
        
        % Load current image
        im = rgb2gray(getFrame(vid_path,v,cFrame));
        
        % Get roi image
        [im_roi,imMask] = isolate_roi(im,x,y,r,phi);
        
        % Stablize rotation
        im_roi = imwarp(im_roi,tform,'OutputView',imref2d(size(im_roi)));
        
        % White out beyond roi
        im_roi(~imMask) = 255;
        
        % Subtract background
        warning off
        im_seq(:,:,i) = uint8(imadjust(imcomplement(imsubtract(im_roi,imRoiMean))));
        warning on
        
        % Visualize
        if 0
            imshow(im_roi2,'InitialMagnification','fit');
            title([num2str(cFrame)])
            pause(0.001);
        end
        
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

if isempty(dir([data_path filesep 'blobParam.mat']))
    
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

for i = 1:length(Rotation)
    
    cIm = im_seq(:,:,i);
    
    % Find blobs in roi
    [props,bw_roi,areas,xB,yB] = findBlobs(cIm,blobParam.tVal,...
                        'area',blobParam.areaMin,blobParam.areaMax);
       
    % Transformation matrix for current roi
    tform = defineSystem2d('roi tform',Rotation(i).roi_rect,...
                          Rotation(i).tform_roi);
    
    % Binary image in global FOR                  
    bw_G = transCoord2d('im L2G',tform,bw_roi,cIm);
    
    if 1
        % Make a truecolor all-green image, make non-blobs invisible
        green = cat(3, zeros(size(im)), ones(size(im)), zeros(size(im)));
        h = imshow(green,'InitialMag','fit');
        set(h, 'AlphaData', bw)
        imshow(im_seq(:,:,i),'InitialMag','fit');
    end
        
end
%TODO: use blobs as masks that are overlaid back onto source image
%TODO: reject masks that move in global FOR




function [im,roi_mask] = isolate_roi(im,x,y,r,theta)

% Maximum size of an image dimension
maxSize = 250;

% rectangular ROI vector
roi_rect = ceil([x-r y-r r*2 r*2]);

% Crop image
im  = imcrop(im,roi_rect);

% Circular coordinates for new roi
xC    = r.*cos(theta)+ size(im,2)/2;
yC    = r.*sin(theta) + size(im,1)/2;

% White out pixels outside of circular roi
roi_mask = roipoly(size(im,1),size(im,2),round(yC),round(xC));

% White out area around roi
im(~roi_mask) = 255;

% Reduce size (helps speed up image registration)
if length(im)>maxSize
    
    imFactor = maxSize/length(im);
    
    im = imresize(im,imFactor);
    
    % White out pixels outside of circular roi
    roi_mask = roipoly(size(im,1),size(im,2),round(yC*imFactor),round(xC*imFactor));
end

% Check for square
if size(im,1)~=size(im,2)
    error('Image not square');
end


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

