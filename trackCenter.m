function varargout = trackCenter(vid_path,v,currDataPath,method,varargin)
% Tracks the motion of an object in a video sequence
% tracker(vid_path,v,method)
%   vid_path - path to video file or image sequence
%   v - structure of info about video (generated by defineVidObject)
%   imInvert - logical that indicates whether to invert the images
%   method - type of tracking approach ('threshold translation' or 'body
%   rotation')
%
% Centroid = tracker(vid_path,v,imInvert,'threshold translation',roi0,tVal,frames,imMean,xMask,yMask) 
% returns just centroid cooridinates
%   roi0 - Initial region-of-interest structure (generated by giveROI('define'))
%   tVal - threshold value
%   frames - listing of frame numbers to analyze (default: all of them)
%   xMask,yMask = coordinates that define a mask 
%
% Centroid = tracker(vid_path,v,imInvert,'thresh trans advanced',iC,frames,imMean,xMask,yMask) 
% uses advanced algoriothm to return just centroid cooridinates
%   iC - Initial conditions structure (with fields for r,x,y,tVal)
%
% Rotation = tracker(vid_path,v,imInvert,'body rotation',roi0,Centroid,frames,imMean) 
% uses advanced algoriothm to return just centroid cooridinates
%   Rotation - Structure with rotation data
%
% Rotation = tracker(vid_path,v,imInvert,'advanced rotation',roi0,Centroid,frames,imMean,xMask,yMask) 
% uses advanced algoriothm to return just centroid cooridinates
%   Rotation - Structure with rotation data

%
%  imInvert - choose 0 for light on dark field, 1 for dark on light field
%  visTracking - Logical to visualize the tracking
%  frames - listing of frame numbers to analyze
%  Centroid - structure with fields x_pix, y_pix, & r_pix for roi coord & radius
%   
%
% Developed by McHenryLab at UC Irvine


%% General parameters

% Number of analyzed frames at which data is save
saveInterval = 10;

% Maximum size of an image dimension (for downsampling)
maxSize = 250;

currDataFile = 'Body.mat';

% Predator & prey colors
clrs{2} = [41 171 226]./255;
clrs{1} = [241 90 36]./255;


%% Load/create data files

if isempty(dir([currDataPath filesep currDataFile]))
    
    error([currDataFile ' (generated by initializeTracking) not found in ' ...
        currDataPath])
else
    % Load Body
    load([currDataPath filesep currDataFile])
end


if isempty(dir([currDataPath filesep 'Initial conditions.mat']))
    error(['Initial conditions.mat (generated by initializeTracking) not found in ' ...
        currDataPath])
else
    % Load iC
    load([currDataPath filesep 'Initial conditions'])
end

% Transfer invert
iminvert = iC.invert;

% Transfer mask
xMask = iC.xTank;
yMask = iC.yTank;



% Interval btwn analyzed frames
anaInterval = iC.frInterval;


%% Parse inputs

if strcmp(method,'centers') 

        imMean = varargin{1};
    

    
%     if length(varargin)>2
%         x0 = varargin{3};
%         y0 = varargin{4};
%         r = varargin{5};
%     else
%         x0 = iC.x;
%         y0 = iC.y;
%         r = iC.r;
%     end

    if length(varargin)>1
        iMode = varargin{2};
    else
        % Interval btwn analyzed frames
        iMode = 1;
    end
    
    if length(varargin)>2
        visSteps = varargin{3};
    else
        visSteps = 0;
    end    
    
    specialAction = 'none';
    
    firstInterval = anaInterval; 
    
    % Number of points to define roi
    numroipts = 400;
    
    propDiff = 0.02;
    
    %visSteps = 1;
    
    clrMode = 'gray';
    %clrMode = 'rgb';
    
    r = iC.r(iMode);
    
    propDiff = iC.propDiff(iMode);
    
    clr = clrs{iMode};
    
    minArea = iC.minArea(iMode);
    maxArea = iC.maxArea(iMode);

elseif strcmp(method,'rotation') 

    % Mean image
    if length(varargin)<1 || iC.useMean==0
        imMean = [];
    else
        imMean = varargin{1};
    end
    
    if length(varargin)>1
        anaInterval = varargin{2};
    else
        % Interval btwn analyzed frames
        anaInterval = 1;
    end
    
    if length(varargin)>5
        visSteps = varargin{6};
    else
        visSteps = 0;
    end    
    
    specialAction = 'none';
    
    firstInterval = anaInterval; 
    
    % Number of points to define roi
    numroipts = 400;
    
    propDiff = iC.propDiff(iMode);
    
    clrMode = 'gray';
    
    
elseif strcmp(method,'refined') 

    % Mean image
    if length(varargin)<1 || ~iC.useMean
        imMean = [];
    else
        imMean = varargin{1};
    end
    
    % Course body data
    courseBod = varargin{2};
    
    r = varargin{3};
    
    if length(varargin)>3
        anaInterval = varargin{4};
    else
        % Interval btwn analyzed frames
        anaInterval = 1;
    end
    
    if length(varargin)>4
        specialAction = varargin{5};
    else
        specialAction = 'none';
    end
    
    if length(varargin)>5
        visSteps = varargin{6};
    else
        visSteps = 0;
    end    

    % Number of points to define roi
    numroipts = 400;
    
    firstInterval = 0;
    
    x0 = courseBod.x(1);
    y0 = courseBod.y(1);
    %idx = 1;

elseif strcmp(method,'arms') 

    % Mean image
    if length(varargin)<1 || ~iC.useMean
        imMean = [];
    else
        imMean = varargin{1};
    end
    
    % Course body data
    courseBod = varargin{2};


    
    if length(varargin)>3
        specialAction = varargin{4};
    else
        specialAction = 'none';
    end
    
    if length(varargin)>4
        visSteps = varargin{5};
    else
        visSteps = 0;
    end    

    % Number of points to define roi
    numroipts = 400;
    
    firstInterval = 0;
    
    x0 = courseBod.x(1);
    y0 = courseBod.y(1);
    
    rMouth = iC.rMouth;
    r = iC.r;


% If no match on method
else
    error('requested method not recognized');
    
end

% Convert mean image to gray for grayscale mode
if strcmp(clrMode,'gray') && size(imMean,3)==3;
    imMean = rgb2gray(imMean);
end


%% Determine frames to analyze
  
% If there are manually-tracked segments . . .
if length(iC.man)>=iMode && ~isnan(iC.man(iMode).start(1))
    
    % Start index
    iStart = 1;
    
    % Set default logic
    ana.auto = logical(zeros(length(Body.frames),1));
    
    % Loop thru manual segments
    for i = 1:length(iC.man(iMode).start)
        
        % Index of end frame
        iEnd = find(Body.frames==iC.man(iMode).start(i),1,'first') - anaInterval;
        
        % Indicies of segment
        idx = iStart:anaInterval:iEnd;
        
        % Designate segment to analyze
        ana.auto(idx) = logical(1);
        
        % New start index
        iStart = iC.man(iMode).end(i) + anaInterval;
    end
    
    % Tack on end of automated portion
    iEnd = length(Body.frames);
    idx = iStart:anaInterval:iEnd;
    ana.auto(idx) = 1;
    
    clear idx iStart iEnd
else
    
    % Analyze all possible frames
    ana.auto = ones(length(Body.frames),1);
    
end

% Figure where to start for center tracking
if strcmp(method,'centers')

    % Indices of last non-nan included in analysis
    idx = find(~isnan(Body.x(:,iMode)) & ana.auto,1,'last');
    
    % Current frame
    cFrame = Body.frames(idx);

    % Starting position 
    x0     = Body.x(idx,iMode);
    y0     = Body.y(idx,iMode);
    %ang0   = Body.ang(idx);
    
end

% Downsmampling used only for image registartion
dSample = 0;

clear idx


%% Smooth data

% if strcmp(method,'refined') || strcmp(method,'arms') 
%     iFrames = ~isnan(courseBod.x);
%     cBod.frames  = courseBod.frames(iFrames);
%     cBod.x       = smooth(courseBod.x(iFrames));
%     cBod.y       = smooth(courseBod.y(iFrames));
%     cBod.ang = smooth(courseBod.ang(iFrames));
%     
%     clear courseBod iFrames
% end


%% Parameter defaults
    
if strcmp(method,'rotation')
    % Initialize image registration parameters
    [optimizer, metric]  = imregconfig('monomodal');
    optimizer.MaximumStepLength = 5e-4;
    optimizer.MaximumIterations = 1500;
    optimizer.RelaxationFactor  = 0.2;
    % optimizer.MaximumStepLength = 10e-4;
    % optimizer.MaximumIterations = 5000;
    % optimizer.RelaxationFactor  = 0.8;
end

% First image
im0 = getFrame(vid_path,v,cFrame,iminvert,clrMode);

%x0 = ana(currSeg).

% Apply mask
im0 = applyMask(im0,xMask,yMask);

% Current roi
roi = giveROI('define','circular',numroipts,r,x0,y0);

% Log first roi
roi0 = roi;

% Focus on roi
%[im0,roi_mask,roi_rect] = isolate_roi(im0,Centroid.x_pix(1),Centroid.y_pix(1),Centroid.r_pix,theta);
[im_roi0,bw_mask] = giveROI('unstabilized',im0,roi,dSample);

% Counter for saving data
nSave = 1;

if strcmp(method,'arms')
    
    % Define roi from approximation
    roi = giveROI('define','circular',numroipts,rMouth,x0,y0);
    
    % Crop translated reference image
    im_roi0 = giveROI('unstabilized',im0,roi,dSample);
    
elseif strcmp(method,'centers')
    
    iFrame = find(Body.frames==cFrame,1,'first');
    
end

% Find blob for starting conditions
if strcmp(method,'centers') 
        
    % If no blob already defined . . .
    if ~isfield(Body,'props') || ...
       (size(Body.props,1)<iFrame) || ...
       size(Body.props,2)<iMode || ...
       isnan(Body.x(iFrame,iMode))
   
        % Find blob at cX,cY
        [props,bwOut] = findBlobs(im0,imMean,propDiff,'coord advanced',x0,y0);
        
        if isstruct(props)
            % Log starting blob
            Body.props(iFrame,iMode).Centroid    = props.Centroid;
            Body.props(iFrame,iMode).PixelList   = props.PixelList;
            Body.props(iFrame,iMode).Area        = props.Area;
            Body.x(iFrame,iMode)                 = props.Centroid(1);
            Body.y(iFrame,iMode)                 = props.Centroid(2);
        else
            error('Blob not found on first frame')
        end
    end
end

% Define frames remaining to analyze
if strcmp(method,'centers')
    % Define vector of frames to analyze
    anaFrames = Body.frames(ana.auto);
    
    % Trim frames already analyzed
    anaFrames = anaFrames(anaFrames>cFrame);
    
    % Starting frame index
    idx = 1;
    
end


%% Tracking object

% Loop thru frames
while true

    
    %TODO: Fix how frames advance
    
    % Index for current frame, one beyond last non-nan
    %idx = find(~isnan(Body.x),1,'last') + anaInterval;
    
    % Stop, if beyond duration
    if idx > length(anaFrames)
        break
    end
    
    % Current frame
    cFrame = anaFrames(idx);   
    
    % Index for current frame
    iFrame = find(Body.frames==cFrame,1,'first');
      
    % Current image
    im = getFrame(vid_path,v,cFrame,iminvert,clrMode);
  
    % Apply mask
    im = applyMask(im,xMask,yMask);
       
    % Index of last non-nan frame
    iLast = find(~isnan(Body.x(:,iMode)) & ana.auto,1,'last');   

    % CENTER TRACKING ------------------------
    if strcmp(method,'centers')
        
        % For predator tracking (uses area filter). . .
        if iMode==2
            minArea = 1;
            maxArea = pi*r^2/4;
        end
        
        % Find blob 
        [props,bwOut] = findBlobs(im,imMean,propDiff,'advanced comparison',...
            Body.props(1,iMode),Body.x(iLast,iMode),Body.y(iLast,iMode),...
            specialAction,minArea,maxArea);
        
%         % For prey tracking (does not use area filter). . .
%         else
%             [props,bwOut] = findBlobs(im,imMean,propDiff,'coord advanced',...
%                 Body.x(iLast,iMode),Body.y(iLast,iMode));
%         end
        
        if isstruct(props)
            %Body.props(iFrame).Centroid    = props.Centroid;
            Body.props(iFrame,iMode).PixelList   = props.PixelList;
            Body.props(iFrame,iMode).Area        = props.Area;
            Body.x(iFrame,iMode)           = props.Centroid(1);
            Body.y(iFrame,iMode)           = props.Centroid(2);
        else
            %Body.props(iFrame).Centroid  = [nan nan];
            Body.props(iFrame,iMode).PixelList = [nan nan];
            Body.props(iFrame,iMode).Area      = nan;
        end
        
        clear props
        
    % REFINED TRACKING ------------------------
    elseif strcmp(method,'refined')
        
        if cFrame < cBod.frames(1)
            currX    = cBod.x(1);
            currY    = cBod.y(1);
            currAng  = cBod.ang(1);
            
        elseif cFrame > cBod.frames(end)
            currX    = cBod.x(end);
            currY    = cBod.y(end);
            currAng  = cBod.ang(end);
            
        else
            
            % Get current approximation
            currX     = interp1(cBod.frames,cBod.x,cFrame);
            currY     = interp1(cBod.frames,cBod.y,cFrame);
            currAng   = interp1(cBod.frames,cBod.ang,cFrame);
        end
        
         % Define roi from approximation
        roi = giveROI('define','circular',numroipts,r,currX,currY);
        
        % Give current image within roi, unstabilized (i.e. unrotated)
        [im_roi,bw_mask_roi] = giveROI('unstabilized',im,roi,dSample);
        
%         % Translate reference image
%         imT = imtranslate(im0,[currX-x0 currY-y0]);
%         
%         % Crop translated reference image
%         im_roicurr = giveROI('stabilized',imT,roi,dSample,-currAng);
  
        % Get properties of blobs
        
        
        tVal = graythresh(im_roi);
        
%         % Get blob
%         bw_roi = im2bw(im_roi,graythresh(im_roi));
%         
%         % Fill holes
%         bw_roi = imfill(~bw_roi,'holes');

        if idx==1
            % Find blob at cX,cY
            [props,bwOut] = findBlobs(im_roi,tVal,'coord advanced',...
                size(im_roi,1)/2,size(im_roi,2)/2,[],specialAction);
            
        else
            % Find blob at cX,cY
            [props,bwOut] = findBlobs(im_roi,tVal,'advanced comparison',...
                size(im_roi,1)/2,size(im_roi,2)/2,Body.props(iFrame-anaInterval),...
                specialAction);
        end
        
        %props = regionprops(bw_roi,'Area','Centroid','Eccentricity','Perimeter');
        
        % Check for first frame
        if idx==1 && length(props)>1         
            error('More than one blob detected on first frame!')
            
        elseif idx==1 && isempty(props)           
            error('No blobs detected on first frame')           
        end           
        
        % If no blobs . . .
        if isempty(props)
            
            warning(['    Frame ' num2str(Body.frames(iFrame)) ': Lost the body!!'])
            
            % Angular rotation up to this point
          %  Body.ang(iFrame)      = nan;
            Body.x(iFrame)        = nan;
            Body.y(iFrame)        = nan;
            Body.props(iFrame)    = nan;
            
        % If blobs . . .
        else     

            % Angular rotation up to this point
          %  Body.ang(iFrame)       = currAng + Dang;
            Body.x(iFrame)         = currX + props.Centroid(1)-size(im_roi,1)/2;
            Body.y(iFrame)         = currY + props.Centroid(2)-size(im_roi,2)/2;
            Body.props(iFrame)     = props;
 
        end
        
        % On first frame .  .
        if idx==1
            x0 = Body.x(1);
            y0 = Body.y(1);
            
            Body.ang(1) = 0;
            
        % In following frames . . .
        else
            
            % Translate reference image
            imT = imtranslate(im0,[Body.x(iFrame)-x0 Body.y(iFrame)-y0]);
            
            % Crop translated reference image
            im_roicurr = giveROI('stabilized',imT,roi,dSample,...
                                            -Body.ang(iFrame-anaInterval));
            
            % Downsample, if necessary
            if length(im_roi)>maxSize
                % Factor by which to resize
                imFactor = maxSize/length(im_roi);
                
                % Downsample images
                im_roiD      = imresize(im_roi,imFactor);
                im_roicurrD = imresize(im_roicurr,imFactor);                 
                
            else
                 % No downsampling
                im_roiD      = im_roi;
                im_roicurrD = im_roicurr;               
            end
            
            % Transformation object to stablize head wrt im0
            tform = imregtform(im_roicurrD,im_roiD,'rigid',optimizer,metric);
            
            % Angular displacement
            Dang   = atan2(tform.T(2,1),tform.T(1,1))*180/pi;
            
            % Log displacement
            Body.ang(iFrame)  = Body.ang(iFrame-anaInterval,1) + Dang;
        end   
   
     % ARM TRACKING ------------------------   
     elseif strcmp(method,'arms')
        
        % Interpolate for current approximate body position
        if cFrame > cBod.frames(end)
            currX    = cBod.x(end);
            currY    = cBod.y(end);
            currAng  = cBod.ang(end);            
        else
            currX     = interp1(cBod.frames,cBod.x,cFrame);
            currY     = interp1(cBod.frames,cBod.y,cFrame);
            currAng   = interp1(cBod.frames,cBod.ang,cFrame);
        end
        
         % Define roi from approximation
        roi = giveROI('define','circular',numroipts,rMouth,currX,currY);
  
        % Give current image within roi, unstabilized (i.e. unrotated)
        [im_roi,bw_mask_roi] = giveROI('unstabilized',im,roi,dSample);
        
        % If first frame, log image
        if idx==1
            
            
            % Store initial values
            Body.x(1)   = x0;
            Body.y(1)   = y0;
            Body.ang(1) = 0;
            Body.tform{1} = affine2d(eye(3));

        else
            
            % Translate reference image to approximate position
            imT = imtranslate(im0,[currX-x0 currY-y0]);
            
            % Crop and rotate translated reference image
            im_roicurr = giveROI('stabilized',imT,roi,dSample,-currAng);
            
            % Downsample, if necessary
            if length(im_roi)>maxSize
                % Factor by which to resize
                imFactor = maxSize/length(im_roi);
                
                % Downsample images
                im_roiD      = imresize(im_roi,imFactor);
                im_roicurrD = imresize(im_roicurr,imFactor);                 
                
            else
                 % No downsampling
                im_roiD     = im_roi;
                im_roicurrD = im_roicurr;               
            end
            
            % Transformation object to stablize head wrt im0
            tform = imregtform(im_roiD,im_roicurrD,'rigid',optimizer,metric,...
                      'InitialTransformation',Body.tform{iFrame-anaInterval});
            
            % Realized scaling factor
            actualFactor = size(im_roiD,1)./ size(im_roi,1);
        
            % Angular & translational displacement
            Dang   = atan2(tform.T(2,1),tform.T(1,1))*180/pi;
            Dx     = tform.T(3,1)./actualFactor;
            Dy     = tform.T(3,2)./actualFactor;
            
            % Store values
            Body.x(iFrame)      = currX - Dx;
            Body.y(iFrame)      = currY - Dy;
            Body.ang(iFrame)    = currAng - Dang;
            Body.tform{iFrame}  = tform;

            clear imT tform 
        end
        
        if visSteps
            im_roi = im_roi0;
            r = rMouth;
        end
    end
    
%     Body.tform{iFrame}      = tform1;

    % Update status
    disp(['bodyTracker (' method ') : done frame ' num2str(Body.frames(iFrame))])   
    
    % Visualize rotation, for debugging 
    if visSteps
        
        % Title text
        t_txt = ['Frame ' num2str(cFrame)];
        
        if strcmp(method,'centers') 
            
            % Start with black image
            bw2 = ~bwOut;
            
            props2 = regionprops(~bwOut,'PixelList');
            
            if strcmp(clrMode,'rgb')
                im2 = applyMask(im,Body.props(iFrame,iMode).PixelList(:,1),...
                    Body.props(iFrame,iMode).PixelList(:,2),[255 255 255],1);
                
                imshow(im2,'InitialMag','fit')
            else
                imshow(im,'InitialMag','fit')
                hold on
                % Make a truecolor all-green image, make non-blobs invisible
                green = cat(3, clr(1).*ones(size(im)), ...
                               clr(2).*ones(size(im)), ...
                               clr(3).*ones(size(im)));
                h = imshow(green,'InitialMag','fit');
                %brighten(bLevel)
                set(h, 'AlphaData', ~bw2)
                
                plot(Body.x(iFrame,iMode),Body.y(iFrame,iMode),'w+')
            end
            
            title(t_txt)
            
            
        else
            
            % Use roi from current frame
            roiCurr = giveROI('define','circular',numroipts,r,Body.x(iFrame),Body.y(iFrame));
            
            %imStable = giveROI('stabilized',im,roi,dSample,tform);
            imStable =  giveROI('stabilized',im,roiCurr,dSample,Body.ang(iFrame));
            theta = linspace(0,2*pi,400);
            
            includeRot = 1;
            
            
            
            % If rotation data included . . .
            
            warning off
            
            subplot(2,2,1)
            imshow(im,'InitialMag','fit')
            % brighten(-0.7)
            hold on
            %plot(Centroid.x(i),Centroid.y(i),'g+',xMask,yMask,'k-')
            line(roiCurr.xPerimG,roiCurr.yPerimG,'Color',[1 0 0 0.2],'LineWidth',3);
            hold off
            title(t_txt)
            
            subplot(2,2,2)
            imshow(imStable,'InitialMag','fit')
            hold on
            %plot(roi.xPerimL,roi.yPerimL,'k-')
            line(roi.xPerimL,roi.yPerimL,'Color',[1 0 0 0.2],'LineWidth',3);
            hold off
            brighten(-0.7)
            
            subplot(2,2,3)
            if idx>1
                imshowpair(im_roi,im_roicurr)
                title('Image comparison')
            end
            
            warning on
        end
        
        
        
        % Pause briefly to render
        pause(0.001)

    end
        
    % Visualize centroid, for debugging 
    if 0
        
        % Title text
        t_txt = ['Frame ' num2str(cFrame) '/' num2str(frames(end))];

        imshow(im,'InitialMag','fit')
        if 1
            brighten(-0.8)
        end
        % brighten(-0.7)
        hold on
        plot(cX,cY,'g+')
        title(t_txt)
        hold off
        
        % Pause briefly to render
        pause(0.001)
    end

     % Clear for next iteration
     clear im_roi tform_roi imStable xC yC h imMask im_roi t_txt

     if nSave > saveInterval
         
         save([currDataPath filesep currDataFile],'Body')
         nSave = 1;
     else
         nSave = nSave + 1;
         
     end
     
     % Advance index
     idx = idx + 1;
     
     clear bwOut bw_mask_roi currX currY currAng cFrame roi
end

% Save data
save([currDataPath filesep currDataFile],'Body')


%% Define outputs

% Threshold method
if strcmp(method,'threshold translation') || ...
        strcmp(method,'threshold roi') || ...
        strcmp(method,'thresh trans advanced') || ...
        strcmp(method,'thresh trans advanced with mask')
    
        varargout{1} = Centroid;
            
% Body rotation method    
elseif strcmp(method,'body rotation') || ...
       strcmp(method,'advanced rotation') || ...
       strcmp(method,'advanced rotation with mask') || ...
       strcmp(method,'advanced rotation with simple mask')
    
    varargout{1} = Rotation;
    
elseif strcmp(method,'visualize') 
    
    if makeVid
        varargout{1} = M;
    else
        varargout{1} = fig;
    end
end



% 
% function im = addmask(im)
% % Mask with distance map (trims fins)
% 
% % Enhance contrast
% im1 = imadjust(im);
% 
% % Binary image
% bw = im2bw(im1,graythresh(im1));
% 
% % Distance map
% bwD = bwdist(bw);
% 
% % Max distance
% maxDist = max(bwD(:));
% 
% % Threshold distance
% threshDist = maxDist/3;
% 
% % Refine binary as above-threshold images
% bw = bwD>threshDist;
% 
% % Dilate, white out outside
% se = strel('disk',3,4);
% bw = imdilate(bw,se);
% im(~bw) = 255;

function im = addsimplemask(im,tVal)
% Adds a regular mask to image
%im1 = imadjust(im);
im1 = im;
bw = im2bw(im1,tVal);
se = strel('disk',3,4);
bw = imdilate(~bw,se);
im(~bw) = 255;

