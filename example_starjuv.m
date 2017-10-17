function example_starjuv

%% Parameters

% Whether to make a video of the data
makeVid = 1;

% Re-run whole analysis
forceReRun = 0;

% Frame rate (fps)
frameRate = 30;

% Calibration constant (m/pix)
calConst = 10^-6; %TODO: Measure this for real

% Execute animation code
do_animate = 0;

% Invert image
imInvert = 1;

% Visualize steps in analysis
visSteps = 0;

% Number of points to define the region of interest
numroipts = 200;


%% Preliminaries

% Extract root directories
paths = givePaths;

% Particular sequence being analyzed
path_seq = ['CSULB juvenile' filesep 'SS18' filesep 'timelapse7'];

% Sample for analysis
vid_path = [paths.vid_root filesep 'Seastars' filesep path_seq];

% Path for data
data_path = [paths.data_root filesep path_seq];

% Load video info (v)
v = defineVidObject(vid_path,'JPG');

% Make data directory, if necessary
if isempty(dir(data_path))
    mkdir(data_path);
end


%% Find mean image

if forceReRun || isempty(dir([data_path filesep 'meanImageData.mat'])) 
    
     % Make mean images 
    imMean = motionImage(vid_path,v,'mean','none',imInvert);
     
    % Save mean image data
    save([data_path filesep 'meanImageData'],'imMean');
    
else
    disp('Loading mean image . . .')
    
    % Load imMean
    load([data_path filesep 'meanImageData.mat'])
end

disp(' ');


%% Interactive mode: Select whether to use mean image

if isempty(dir([data_path filesep 'Initial conditions.mat'])) || forceReRun
    
    % First image, with mean image subtraction
    im0Mean = getFrame(vid_path,v,v.UserData.FirstFrame,imInvert,'gray',imMean);
    
    % First image
    im0NoMean = getFrame(vid_path,v,v.UserData.FirstFrame,imInvert,'gray');
    
    f = figure;
    subplot(1,2,1)
    imshow(imadjust(im0Mean),'InitialMag','fit')
    title('Mean image subtracted')
    
    subplot(1,2,2)
    imshow(imadjust(im0NoMean),'InitialMag','fit')
     title('No mean image subtracted')
    
    b = questdlg('Are the tube feet more visible on the right or left image?',...
        '','Left','Right','Cancel','Right');
    
    if strcmp(b,'Cancel')
        return
        
    elseif strcmp(b,'Left')
        
        im = im0Mean;
        useMean = 1;
        
    elseif strcmp(b,'Right')   
        
        imMean = [];
        useMean = 0;
        
        % Overwrite image image data
        save([data_path filesep 'meanImageData'],'imMean');
        
        im = im0NoMean;
    end
    
    close(f);
    
    % Initial position
    disp(' ')
    disp('Select animal to be tracked')
    [x,y] = imInteract(im,'points',1);
    
    % Threshold
    disp(' ')
    disp('Select threshold')
    tVal = imInteract(im,'threshold');
    
    % Radius
    disp(' ')
    disp('Select roi radius')
    r = imInteract(im,'radius',x,y);
    
    % Store data
    iC.x       = x;
    iC.y       = y;
    iC.tVal    = tVal;
    iC.r       = r;
    iC.useMean = useMean;
    
    % Save data
    save([data_path filesep 'Initial conditions'],'iC')
    
    clear im useMean im0Mean im0NoMean x y tVal r
else
    
    disp('Loading initial condition data . . .')

    % Load initial conditions, iC
    load([data_path filesep 'Initial conditions.mat']);
    
    % Load imMean
    load([data_path filesep 'meanImageData.mat'])
    
end

disp(' ')
       

%% Track centroid coordinates

if isempty(dir([data_path filesep 'Centroid.mat'])) || forceReRun
    
    % Region of interest for first frame
    roi0 = giveROI('define','circular',numroipts,iC.r,iC.x,iC.y);
    
    % Run tracker code for centroid
    Centroid = tracker(vid_path,v,imInvert,'threshold translation',...
        roi0,iC.tVal);
    
    % Save data
    save([data_path filesep 'Centroid'],'Centroid')
    
    close
else
    
    disp('Loading Centroid coordinates . . .')

    % Load 'Centroid'
    load([data_path filesep 'Centroid.mat']);
    
end

disp(' ')


%% Track rotation

if isempty(dir([data_path filesep 'Rotation.mat'])) || forceReRun
    
    % Region of interest for first frame
    roi0 = giveROI('define','circular',numroipts,iC.r,iC.x,iC.y);
    
    % Run tracker code for centroid
    Rotation = tracker(vid_path,v,imInvert,'body rotation',roi0,Centroid,imMean);
    
    % Save data
    save([data_path filesep 'Rotation'],'Rotation')
    
else
    % Load 'Rotation'
    load([data_path filesep 'Rotation.mat']);
end


%% Store coordinate transformations in S structure

% Region of interest for first frame
roi0 = giveROI('define','circular',numroipts,iC.r,iC.x,iC.y);

% Create coordinate transformation structure 
S = defineSystem2d('roi',roi0,Centroid,Rotation);

clear iC Centroid Rotation roi0


%% Visualize and approve results up to this point

if 0
    % Make movie of tracking up to this point
    M = aniData(vid_path,v,imInvert,'Centroid & Rotation',S,1);
    
    while true
        
        % Prompt to proceed
        ButtonName = questdlg('Does the tracking look good?', ...
            '', 'Yes, proceed', 'No, Cancel', 'Replay animation', 'Yes, proceed');
        
        switch ButtonName,
            case 'No, Cancel',
                disp(''); disp('Quitting analysis . . .');beep
                return
                
            case 'Yes, proceed',
                
                clear M 
                
                break
        end % switch
        
        clear ButtonName
        
        movie(M)
    end
end


%% Make mean image of roi 

if forceReRun || isempty(dir([data_path filesep 'meanRoiImageData.mat']))
    
     % Make mean images 
    imRoiMean = motionImage(vid_path,v,'mean roi','none',imInvert,100,S,0);
    
    % Plot mean images
    f = figure;

    imshow(imRoiMean,'InitialMag','fit');
    title('Mean roi image')
    
    % Prompt to proceed
    ButtonName = questdlg('Does the roi mean image look good?', ...
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
    
else
    disp('Loading mean images . . .')
    
    % Load imRoiMean
    load([data_path filesep 'meanRoiImageData.mat'])
    
end

disp(' ');


%% Interactive mode: select threshold and area limits

if forceReRun || isempty(dir([data_path filesep 'blobParam.mat']))
    
    justReplay = 0;
    
    f = figure;
    
    dSample = 0;

    while true
        
        if ~justReplay
            
            i = 1;
            
            % Full-frame image
            im = getFrame(vid_path,v,S.frames(i),imInvert,'gray');  
            
            % Roi image, mean image subtracted
            im_roi = giveROI('stabilized',im,S.roi(i),dSample, ...
                S.tform(:,:,i),imRoiMean);
            
            % Find threshold
            tVal = imInteract(im_roi,'threshold');
            %tVal = 0.8349;
            
            % Find area
            [areaMin,areaMax] = imInteract(im_roi,'area',tVal);
            %areaMin = 996.4110;areaMax = 2.3445e+04;
            
            M = aniData(vid_path,v,imInvert,'blobs L simple',S,tVal,areaMin,...
            areaMax,1,imRoiMean,dSample);
        
        else
            movie(M)
        end
        
        % Reset
        justReplay = 0;
        
        % Prompt
        choice = questdlg('Do threshold and area values look good?', ...
            '', 'Yes, proceed','No, redo','Replay','Yes, proceed');
        
        % Break, if good
        if strcmp(choice,'Yes, proceed')
            
            % Store parameters
            blobParam.tVal    = tVal;
            blobParam.areaMin = areaMin;
            blobParam.areaMax = areaMax;
            blobParam.AR_max  = 1.5;
            
            % Close figure window
            close(f)
            
            clear M tVal areaMin areaMax im im_roi i
            
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

if forceReRun || isempty(dir([data_path filesep 'Blob data.mat']))
    
    % Downsample
    dSample = 0;
    
    %blobParam.AR_max  = 1.5;
    
    % Get blobs
    B = anaBlobs(vid_path,v,'G&L props',S,blobParam,imInvert,dSample,...
                 imRoiMean);
    
    % Save image image data
    save([data_path filesep 'Blob data'],'B');
    
else
    
    disp('Loading B structure . . .');
    
    % Load 'B' structure
    load([data_path filesep 'Blob data']);
end

disp(' ');


%% Interactive mode: Select frame interval

if forceReRun || isempty(dir([data_path filesep 'winDur.mat']))
    
    % Whether to downsample the image
    dSample = 0;
    
    f = figure;
    
    winDur = 10;
    
    % Frame index
    idx = 2+ceil(winDur/2);
    
    while true
        
        % Starting and ending frames
        frStart = S.frames(idx-floor(winDur/2));
        frEnd   = S.frames(idx+ceil(winDur/2));
        
        % Image of static motion
        imB = motionImage(vid_path,v,'bw static',[frStart:frEnd],B);
        
        % Full-frame image
        im1 = getFrame(vid_path,v,frStart,imInvert,'gray');
        
        % Full-frame image
        im2 = getFrame(vid_path,v,frEnd,imInvert,'gray');
        
        subplot(1,2,1)
        imshowpair(im1,im2)
        title(['Window length = ' num2str(winDur)])
        
        subplot(1,2,2)
        imshow(imB,'InitialMag','fit')
        title(['Window length = ' num2str(winDur)])
 
        % Prompt
        b = questdlg('Does the window length look good?', '', ...
                         'Yes, proceed', 'No, set', 'Cancel', 'Yes, proceed');
        if strcmp(b,'Yes, proceed')
            close(f)
            break
            
        elseif strcmp(b,'No, set')
            an = inputdlg({'Window length (frames)'},'',1,{num2str(winDur)});
            winDur = str2num(an{1});

        else
            return
        end

        clear im1 im2 frStart frEnd b an 
    end
    
    % Prompt for threshold
    winThresh = imInteract(imB,'threshold');
    
    % Save winDur
    save([data_path filesep 'winDur'],'winDur');  
    
    % Save winDur
    save([data_path filesep 'winThresh'],'winThresh');  
    
else
    disp('Loading winDur . . .');
    
    % Load winDur, parameter values
    load([data_path filesep 'winDur.mat'])
    
    % Load winThresh, parameter values
    load([data_path filesep 'winThresh.mat'])
end

disp(' ');


%% Filter out moving blobs

if 1 %forceReRun || isempty(dir([data_path filesep 'Blob, filtered data.mat']))

    % Find blobs that don't move in global FOR
    Bf = anaBlobs(vid_path,v,'filter motion',B,winDur,S,blobParam,winThresh);
    
    % Save blob data
    save([data_path filesep 'Blob, filtered data'],'Bf','-v7.3');
    
else
     disp('Loading Bf structure . . .');
    
    % Load 'Bf' structure
    load([data_path filesep 'Blob, filtered data']);
end

disp(' ');disp(' ')


%% Display results

disp('Creating movie . . .')
M = aniData(vid_path,v,imInvert,'blobs G&L',Bf,0);

rrr=2

%% Isolate tube feet

%isolateRoiMotion(vid_path,v,data_path,Rotation,Centroid,imInvert)


return

%% Calculate stuff

% Frame period
dt = 1/frameRate;

% Step thru data
for i = 1:length(Rotation)
    
    % Time
    d.t(i,1) = i.*dt;
    
    % Head angle correction from image registration (imStable)
    d.theta(i,1) = atan2(Rotation(i).tform_roi.T(1,2),Rotation(i).tform_roi.T(1,1));
    
    % Coordinate data
    d.x(i,1) = Centroid.x_pix(i) * calConst;
    d.y(i,1) = Centroid.y_pix_flip(i) * calConst;
    
end

subplot(2,1,1)
plot(d.t,d.x.*1000,'-',d.t,d.y.*1000,'-')
xlabel('time (s)')
ylabel('Position (mm)')
legend('x','y')

subplot(2,1,2)
plot(d.t,d.theta.*180/pi,'-')
ylabel('Orientation (deg)')


%% Animation 

if do_animate
    
    % Make movie
    M = tracker(vid_path,v,'visualize',blobMultiple,0,Centroid.frames,Centroid,...
        Rotation,makeVid);
    
    % Write movie to disk
    if makeVid
        vid_save_path = uigetdir(data_path,'Save movie');
        vInfo = VideoWriter([vid_save_path filesep 'video.mp4'],'MPEG-4');
        vInfo.FrameRate = 15;
        open(vInfo)
        writeVideo(vInfo,M)
        close(vInfo)
    end
    
end

 return
% 
% % Frame numbers
% frames = v.UserData.FirstFrame:v.UserData.LastFrame;
% 
% % First image
% im = rgb2gray(getFrame(vid_path,v,frames(1)));
% 
% % Select sea star from first frame
% [tVal,x,y] = imInteract(im,'threshold and selection');
% pause(0.5)
% close
% 
% % Angular coordinates for roi
% theta = linspace(0,2*pi,500);
% 
% 
% %% Define reference image and blob size
% 
% % Define binary image of seastar
% bw = im2bw(im,tVal);
% bw = bwselect(bw,x,y);
% bw = imfill(bw,'holes');
% 
% % Find perimeter
% [yPerim, xPerim] = find(bwperim(bw,8));
% 
% % Blob diameter
% blobDiam = max([range(xPerim) range(yPerim)]);
% 
% % Circular coordinates for roi
% r     = blobMultiple*blobDiam/2;
% xC    = r.*cos(theta)+x;
% yC    = r.*sin(theta)+y;
% 
% % Bitmap mask of roi
% bwMask = roipoly(bw.*0,xC,yC);
% 
% % Replicate source image and mask outside
% imRef = im;
% imRef(~bwMask) = 255;
% 
% % Clear unneeded
% clear r xC yC xPerim yPerim bwMask bw im
% 
% 
% %% Track thru frames
% 
% % Initialize image registration
% [optimizer, metric]  = imregconfig('monomodal');
% optimizer.MaximumStepLength = 5e-4;
% optimizer.MaximumIterations = 1500;
% optimizer.RelaxationFactor  = 0.2;
% 
% % Loop thru frames
% for i = 5 %1:length(frames)
%     
%     % Current frame
%     cFrame = frames(i);
%     
%     % Current image
%     im = rgb2gray(getFrame(vid_path,v,cFrame));
% 
%     % Make blob image 
%     bw = im2bw(im,tVal);
%     bw = bwselect(bw,x,y);
%     bw = imfill(bw,'holes');
%     
%     % Find center of blob
%     props = regionprops(bw,'Centroid');
%     
%     % Check for one (and only one) blob        
%     if length(props)>1    
%         error('More than one object selected');
%         
%     elseif isempty(props)
%         error('No object selected');
%     end
%     
%     % Find new centroid
%     x = props.Centroid(1);
%     y = props.Centroid(2);
%     
%     % Circular coordinates for roi 
%     r     = blobMultiple*blobDiam/2;
%     xC    = r.*cos(theta)+x;
%     yC    = r.*sin(theta)+y;
%     
%     % Bitmap mask of roi
%     bwMask = roipoly(bw.*0,xC,yC);
%     
%     % Mask outside
%     im(~bwMask) = 255;
%     
%     % Transformation matrix for current image
%     tformC = imregtform(im,imRef,'rigid',optimizer, metric);
%     
%     imHead = imwarp(im,invert(tformC),'OutputView',imref2d(size(im)));
%     
%     % Clear of next loop
%     clear r xC yC props bw bwMask im cFrame
% end
% 
% return





