function initializeTracking(vid_path,v,currDataPath,method,varargin)
% Creates data files for analysis by kineBox

%TODO: Add pred-prey mode

% Check for path in data dir
if isempty(dir(currDataPath))
    % Make directory, if not there
    mkdir(currDataPath);
end

%% Parse inputs

% Mean image
if length(varargin)>0
    imMean = varargin{1};
else
    imMean = [];
end

% Image invert (default)
if length(varargin)>1
    imInvert = varargin{2};
else
    imInvert = 0;
end

% Image invert (default)
if length(varargin)>1
    imInvert = varargin{2};
else
    imInvert = 0;
end


%% Interactively select initial conditions

if isempty(dir([currDataPath filesep 'Initial conditions.mat']))
    
    % DURATION -------------------------------------
    
    % Get frame duration
    iC = selectDuration(vid_path,v,imMean);
    
    answer = inputdlg({'Interval btwn frames for COURSE analysis',...
                       'Interval btwn frames for FINE analysis'},'',1,{'10','1'});
    
    % Interval between tracked frames
    iC.frIntervalCourse = str2num(answer{1});
    iC.frIntervalFine   = str2num(answer{2});
     
    % First image
    im = getFrame(vid_path,v,iC.startFrame,imInvert,'gray',imMean);
    
    % IMAGE PROPERTIES -------------------------------------
    
    % First image (No mean-image subtraction)
    imNoMean = getFrame(vid_path,v,...
                iC.startFrame,imInvert,'gray');
    
    figure;
    subplot(1,2,1)
    imshow(imNoMean); title('No mean-image subtraction')
    subplot(1,2,2)
    imshow(imadjust(im)); title('With mean-image subtraction')
            
    ButtonName = questdlg('Use mean-image subtraction?', ...
        '','Yes', 'No', 'Cancel', 'Yes');
    
    if strcmp(ButtonName,'Yes')
        iC.useMean = 1;
        
        % Brightness level (darkens images)
        bLevel = -0.5;
        
    elseif strcmp(ButtonName,'No')
        iC.useMean = 0;
        
        imMean = [];
        
        bLevel = 0;
        
    else
        return
    end
    
    close
    clear buttonName 
    
    
    imI = getFrame(vid_path,v,iC.startFrame,0,'gray',imMean);
    
    figure;
    imshow(imI,'InitialMag','fit')
    
    ButtonName = questdlg('Is the body dark or light?', ...
        '','Dark', 'Light', 'Cancel', 'Dark');
    
    if strcmp(ButtonName,'Dark')
        iC.invert = 0;
        
    elseif strcmp(ButtonName,'Light')
        iC.invert = 1;
        
    else
        return
    end
    
    clear buttonName imI
    close
    
    
       % First image (No mean-image subtraction)
    imNoMean = getFrame(vid_path,v,...
                iC.startFrame,iC.invert,'gray');
       
    
    % MASKING -------------------------------------
    
    ButtonName = questdlg('Were the experiments run in an arena?', ...
        '','Yes, elliptical', 'Yes, rectangular','No', 'Yes, elliptical');
    
    if strcmp(ButtonName,'Yes, elliptical')
        
        % Boundaries of tank
        disp(' ')
        disp('Select elliptical boundaries of the tank')
        disp(' ')
        [iC.xTank,iC.yTank] = imInteract(imNoMean,'ellipse');
    
    elseif strcmp(ButtonName,'Yes, rectangular')
        
        % Boundaries of tank
        disp(' ')
        disp('Select rectangular boundaries of the tank')
        disp(' ')
        
        %TODO: Make this mode:
        [iC.xTank,iC.yTank] = imInteract(imNoMean,'rectangle');
        
    elseif strcmp(ButtonName,'No')
        iC.xTank = [];
        iC.yTank = [];
        
    else
        return
    end
    
    clear buttonName 
    
    % CALIBRATION -------------------------------------

    figure;
    imshow(imNoMean,'InitialMag','fit')
    
     ButtonName = questdlg('Are there calibration landmarks in view?', ...
        '','Yes', 'No', 'Cancel', 'Yes');
    
    if strcmp(ButtonName,'Yes')
        
        % Calibration
        disp(' ')
        disp('Select two points for spatial calibration')
        [xCal,yCal] = imInteract(imNoMean,'points',2);
        
        answer = inputdlg('Distance between 2 points (m)?','',1,{'0.5'});

        iC.xCalPts   = xCal;
        iC.yCalPts   = yCal;
        iC.calconst = str2num(answer{1}) / hypot(diff(xCal),diff(yCal));
        
    elseif strcmp(ButtonName,'No')
        iC.xCal = [];
        iC.yCal = [];
        
    else
        return
    end
    
    clear buttonName xCal yCal answer
    
    
    % FRAME RATE ----------------------------------------
    
    if isfield(v,'FrameRate')

        iC.frameRate = v.FrameRate;
    else
        answer = inputdlg({'Frame rate (fps)'},'',1,{'29.97'});
        
        iC.frameRate = str2num(answer{1});
    end
    
    clear answer
            
    
    % SELECT ANIMAL POSITION ------------------------------
    
    % First image
    im = getFrame(vid_path,v,iC.startFrame,iC.invert,'gray',imMean);
         
    % Single body tracking -----------
    if strcmp(method,'single body') || strcmp(method,'single body, with arms') 
            
        % Initial position
        disp(' ')
        disp('Select animal to be tracked')
        [iC.x,iC.y] = imInteract(im,'points',1);
        
        % Radius
        disp(' ')
        disp('Select roi radius')
        iC.r = imInteract(im,'radius',iC.x,iC.y);
        
   
    % Two body tracking (i.e. predator-prey) ---
    elseif strcmp(method,'pred prey')
        
        % Initial position of predator
        disp(' ')
        disp('Select predator centroid')
        [iC.xPd0,iC.yPd0] = imInteract(im,'points',1,bLevel);
        
        % Radius for predator
        disp(' ')
        disp('Select roi radius for predator')
        iC.rPd = imInteract(im,'radius',iC.xPd0,iC.yPd0,bLevel);
        
        % Initial position of prey
        disp(' ')
        disp('Select prey centroid')
        [iC.xPy0,iC.yPy0] = imInteract(im,'points',1,bLevel);
        
        % Radius for prey
        disp(' ')
        disp('Select roi radius for prey')
        iC.rPy = imInteract(im,'radius',iC.xPy0,iC.yPy0,bLevel);
        
    end
     
    % Arms
    if strcmp(method,'single body, with arms') 
        % Prompt for arms
        answer = inputdlg({'Number of arms'},'',1,{'5'});
         
        % Log response
        iC.armNum = str2num(answer{1});     
        
        disp(' ')
        disp('Select tip of arms centroid')
        
        % Capture coordinates of arm tips
        [iC.xArms,iC.yArms] = imInteract(im,'points',iC.armNum,bLevel);  
        
        % Radius around mouth
        disp(' ')
        disp('Select radius around body area (without arms)')
        iC.rMouth = imInteract(im,'radius',iC.x,iC.y);
    end
    
    % Save initial conditions
    save([currDataPath filesep 'Initial conditions'],'iC')
    
    clear im
    
else
    load([currDataPath filesep 'Initial conditions'])
end


%% Make data files

if strcmp(method,'single body') && isempty(dir([currDataPath filesep 'Body course.mat']))
      
    % Make empty Centroid structure
    Body.frames     = [iC.startFrame:iC.endFrame]';
    Body.tform{1}   = affine2d(eye(3));
    Body.x          = nan(size(Body.frames));
    Body.y          = nan(size(Body.frames));
    Body.ang        = nan(size(Body.frames));
    %Body.y_flip     = nan(size(Body.frames));
    
    % Initial values
    Body.x(1)        = iC.x;
    Body.y(1)        = iC.y;
    Body.ang(1)      = 0;
    
    % Save data
    save([currDataPath filesep 'Body course'],'Body')
    
elseif isempty(dir([currDataPath filesep 'Body arm.mat']))
      
    % Make empty Centroid structure
    Body.frames     = [iC.startFrame:iC.endFrame]';
    Body.x          = nan(size(Body.frames));
    Body.y          = nan(size(Body.frames));
    Body.ang        = nan(size(Body.frames));
    %Body.y_flip     = nan(size(Body.frames));
    
    % Initial values
    Body.x(1)        = iC.x;
    Body.y(1)        = iC.y;
    Body.ang(1)      = 0;
    
    % Save data
    save([currDataPath filesep 'Body arm'],'Body')
    
    % Two body tracking (i.e. predator-prey) ---
elseif strcmp(method,'pred prey') && ...
        isempty(dir([currDataPath filesep 'Pd Body course.mat']))
    
    % Make empty Centroid structure
    Body.frames     = [iC.startFrame:iC.endFrame]';
    Body.tform{1}   = affine2d(eye(3));
    Body.x          = nan(size(Body.frames));
    Body.y          = nan(size(Body.frames));
    Body.ang        = nan(size(Body.frames));
    %Body.y_flip     = nan(size(Body.frames));
    
    % Initial values
    Body.x(1)        = iC.xPd0;
    Body.y(1)        = iC.yPd0;
    Body.ang(1)      = 0;
    
    % Save data
    save([currDataPath filesep 'Pd Body course'],'Body')
    save([currDataPath filesep 'Pd Body fine'],'Body')
    
    % Initial values
    Body.x(1)        = iC.xPy0;
    Body.y(1)        = iC.yPy0;
    Body.ang(1)      = 0;
    
    save([currDataPath filesep 'Py Body course'],'Body')
    save([currDataPath filesep 'Py Body fine'],'Body')

end





function clipInfo = selectDuration(vidPath,v,imMean)
% Interactively prompts to select a duration for analysis

% Make figure
f = figure;

% Default first and end frames
firstFrame = v.UserData.FirstFrame;
lastFrame  = v.UserData.LastFrame;

% Loop for multiple tries at frame numbers
while true

    % Get images
    im1 = getFrame(vidPath,v,firstFrame,0,'gray',imMean);
    im2 = getFrame(vidPath,v,lastFrame,0,'gray',imMean);
    
    % Plot candidate frames
    subplot(1,2,1)
    imshow(im1,'InitialMag','fit')
    title(['First frame (' num2str(firstFrame) ')'])
    
    subplot(1,2,2)
    imshow(im2,'InitialMag','fit')
    title(['Last frame (' num2str(lastFrame) ')'])
    
    % Ask for approval on frames
    an = questdlg('Are these good starting and ending frames?','','Yes',...
        'No','Cancel','Yes');
    
    if strcmp(an,'Yes')
        
        close(f)
        break
        
    elseif strcmp(an,'No')
        
        prompt={'Start frame num:','Last frame num:'};
        name='Choose clip duration';
        numlines=1;
        defaultanswer={num2str(firstFrame),num2str(lastFrame)};
        
        answer = inputdlg(prompt,name,numlines,defaultanswer);
        
        if isempty(answer)
            return
        end
        
        firstFrame   = str2num(answer{1});
        lastFrame    = str2num(answer{2});
        
    else
        return
    end
end

clipInfo.startFrame = firstFrame;
clipInfo.endFrame   = lastFrame;



