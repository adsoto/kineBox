function initializeTracking(vid_path,v,currDataPath,method,varargin)


%% Parameters

% Currently only support grayscale mode
clrMode = 'gray';

% Predator & prey colors
pyClr = [41 171 226]./255;
pdClr = [241 90 36]./255;

% Level to brighten the movie
bLevel = 0.5;


%% Create data files for analysis by kineBox

% Check for path in data dir
if isempty(dir(currDataPath))
    % Make directory, if not there
    mkdir(currDataPath);
end


%% Preliminary questions

if isempty(dir([currDataPath filesep 'Initial conditions.mat']))
    
     ButtonName = questdlg('Which kind of experiment?', ...
                          '','Single body', 'Pred-prey', 'Cancel', 'Single body');
                      
    if strcmp(ButtonName,'Single body')
        iC.expType = 'single';
        
    elseif strcmp(ButtonName,'Pred-prey')
        iC.expType = 'pred prey';
    else
        return
    end
    
    % Frame rate & interval for analysis
    if isfield(v,'FrameRate')

        iC.frameRate = v.FrameRate;
        
        answer = inputdlg({'Interval btwn frames for analysis (frames)'},'',1,{'1'});
        
        % Interval between tracked frames
        iC.frInterval = str2num(answer{1});
    else
        answer = inputdlg({'Interval btwn frames for analysis (frames)',...
                           'Frame rate (fps)'},'',1,{'1','29.97'});
        
        iC.frInterval  = str2num(answer{1});
        iC.frameRate   = str2num(answer{2});
    end 
    
    % Manual analysis
    bName = questdlg('Any intervals to analyze manually?', ...
                          '','Yes', 'No', 'Cancel', 'Yes');
                      
    if strcmp(bName,'Cancel') || isempty(bName)
        return
        
    elseif strcmp(bName,'Yes')
        
        % Get number of intervals
        answer = inputdlg({'How many intervals?'},'',1,{'1'}); 
        numInt = str2num(answer{1});
        
        % Start with empty data
        iC.man(1).start = [];
        iC.man(1).end   = [];
        
        if strcmp(iC.expType,'pred prey')
            iC.man(2).start = [];
            iC.man(2).end   = [];
        end
        
        % Loop thru intervals
        for i = 1:numInt

            if strcmp(iC.expType,'pred prey')
                answer = inputdlg({'Start frame?','End frame','Pred (P), Prey (E), or both (PE)?'},...
                    ['Interval ' num2str(i)],1,{'','','P'});
                
                % Log predator interval
                if strcmp(answer{3},'P') || strcmp(answer{3},'p')
                    iC.man(1).start  = [iC.man(1).start; str2num(answer{1})];
                    iC.man(1).end    = [iC.man(1).end; str2num(answer{2})];
                    
                    % Log prey interval
                elseif strcmp(answer{3},'E') || strcmp(answer{3},'e')
                    iC.man(2).start  = [iC.man(2).start; str2num(answer{1})];
                    iC.man(2).end    = [iC.man(2).end; str2num(answer{2})];
                    
                    % Log interval for both predator and prey
                elseif strcmp(answer{3},'PE') || strcmp(answer{3},'pe')||...
                        strcmp(answer{3},'EP') ||  strcmp(answer{3},'ep')
                    iC.man(1).start  = [iC.man(1).start; str2num(answer{1})];
                    iC.man(1).end    = [iC.man(1).end; str2num(answer{2})];
                    iC.man(2).start  = [iC.man(2).start; str2num(answer{1})];
                    iC.man(2).end    = [iC.man(2).end; str2num(answer{2})];
                    
                else
                    error('Input not recognized');
                end
                
            else
                answer = inputdlg({'Start frame?','End frame'},...
                    ['Interval ' num2str(i)],1,{'','','P'});
                
                % Log interval
                iC.man(1).start  = [iC.man(1).start; str2num(answer{1})];
                iC.man(1).end    = [iC.man(1).end; str2num(answer{2})];
            end
        end
        
    elseif strcmp(bName,'No')
        
        iC.man(1).start  = nan;
        iC.man(1).end    = nan;
    end
    
    clear bName answer
end


%% Choose movie duration

if isempty(dir([currDataPath filesep 'Initial conditions.mat']))
    
    % Make figure
    f = figure;
    
    % Default first and end frames
    firstFrame = v.UserData.FirstFrame;
    lastFrame  = v.UserData.LastFrame;
    
    % Loop for multiple tries at frame numbers
    while true
        
        % Get images
        im1 = getFrame(vid_path,v,firstFrame,0,'rgb');
        im2 = getFrame(vid_path,v,lastFrame,0,'rgb');
        
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
    
    % Store start and end frames
    iC.startFrame = firstFrame;
    iC.endFrame   = lastFrame;

    clear im1 im2 f prompt name numlines answer firstFrame lastFrame
end


%% Determine whether to invert image

if isempty(dir([currDataPath filesep 'Initial conditions.mat']))
    
    im = getFrame(vid_path,v,iC.startFrame,0,'gray');
    
    f = figure;
    imshow(im,'InitialMag','fit')
    
    ButtonName = questdlg('Is the foreground body dark or light?', ...
                          '','Dark', 'Light', 'Cancel', 'Dark');
    
    if strcmp(ButtonName,'Dark')
        iC.invert = 0;
        
    elseif strcmp(ButtonName,'Light')
        iC.invert = 1;
        
    else
        return
    end
    
    clear buttonName im
    close(f)
    
end


%% Tank mask       
    
if isempty(dir([currDataPath filesep 'Initial conditions.mat']))  

    im = getFrame(vid_path,v,iC.startFrame,iC.invert,'gray');
    
    f = figure;
    imshow(im,'InitialMag','fit')
    
    ButtonName = questdlg('Were the experiments run in an arena?', ...
        '','Yes, elliptical', 'Yes, rectangular','No', 'Yes, elliptical');
    
    if strcmp(ButtonName,'Yes, elliptical')
        
        % Boundaries of tank
        disp(' ')
        disp('Select elliptical boundaries of the tank')
        disp(' ')
        [iC.xTank,iC.yTank] = imInteract(im,'ellipse');
    
    elseif strcmp(ButtonName,'Yes, rectangular')
        
        % Boundaries of tank
        disp(' ')
        disp('Select rectangular boundaries of the tank')
        disp(' ')
        
        %TODO: Make this mode:
        [iC.xTank,iC.yTank] = imInteract(im,'rectangle');
        
    elseif strcmp(ButtonName,'No')
        iC.xTank = [];
        iC.yTank = [];
        
    else
        return
    end
    
    close(f)
    clear buttonName im
end


%% Calculate mean image

if isempty(dir([currDataPath filesep 'meanImageData.mat']))
    
    % Calculate mean image
    imMean = motionImage(vid_path,v,'mean bright','none',iC.invert);
    
    if iC.invert==1
        imMean = imcomplement(imMean);
    end
    
    % Save mean image data
    save([currDataPath filesep 'meanImageData'],'imMean');
    
else
    % Load imMean
    load([currDataPath filesep 'meanImageData.mat'])
end
    

%% Determine contrast difference between blobs and mean image

if isempty(dir([currDataPath filesep 'Initial conditions.mat']))    
    
    % Frames for plot
    visFrames = round(linspace(iC.startFrame,iC.endFrame,9));
    
    % Make figure
    f = figure;
    
    % Render 9 frames with current settings
    for i = 1:9
        % Get image, subtract around tank
        im{i} = getFrame(vid_path,v,visFrames(i),iC.invert,'gray');
        im{i} = applyMask(im{i},iC.xTank,iC.yTank);
        
        subplot(3,3,i)
        imshow(im{i},'InitialMag','fit')
        title(['Frame ' num2str(visFrames(i))])
    end
    
    brighten(bLevel)
    
    % Index for run thru loop
    runNum = 1;
    
    % Loop thru for each body
    while true
     
        % Starting value
        propDiff = 0.15;

        disp(' ')
        
        if strcmp(iC.expType,'pred prey') && runNum==1
            disp('PREDATOR ----------------------------') 
            clr = pdClr;
            
        elseif strcmp(iC.expType,'pred prey')
            disp('PREY --------------------------------')      
            clr = pyClr;
            
        else
            clr = [0 1 0];
        end
        
        disp('Interactively adjust proportional difference and area factor');
        disp(' ')
        disp('Up arrow :      More inclusive in contrast')
        disp('Down arrow :    Less inclusive in contrast')
        disp(' ')
        
        while true
            
            % Adjust blob areas
            for i = 1:9
                % Find blob at cX,cY
                [props,bwOut] = findBlobs(im{i},imMean,propDiff,'all');
                
                figure(f);
                
                subplot(3,3,i)
                hold on
                % Make a truecolor all-green image, make non-blobs invisible
                clrField = cat(3, clr(1).*ones(size(im{i})), ...
                                  clr(2).*ones(size(im{i})), ...
                                  clr(3).*ones(size(im{i})));
                h(i) = imshow(clrField,'InitialMag','fit');
                %brighten(bLevel)
                set(h(i), 'AlphaData', bwOut)
                hold off
            end
            
            % Interactive mode
            [x,y,b] = ginput(1);
            
            % Up arrow
            if b==30
                propDiff = propDiff/1.25;
                
                % Down arrow
            elseif b==31
                propDiff = propDiff*1.25;
                
                
            elseif isempty(b)
                iC.propDiff(runNum,1) = propDiff;
                break
            end
            
            disp(['      propDiff = ' num2str(propDiff)])
            
            delete(h)
        end
        
        % If first of 2 iterations
        if strcmp(iC.expType,'pred prey') && runNum==1
            runNum = 2;
            
        % Otherwise, leave
        else
            close(f)
            clear frames im
            break
        end
    end
    clear b bwOut clr clrField f h i 
end  


%% Animal position(s)

if isempty(dir([currDataPath filesep 'Initial conditions.mat']))   
    
    % First image
    im = getFrame(vid_path,v,iC.startFrame,iC.invert,'rgb');
         
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
        [iC.x,iC.y] = imInteract(im,'points',1);
        
        % Radius for predator
        disp(' ')
        disp('Select roi radius for predator')
        iC.r = imInteract(im,'radius',iC.x,iC.y);
        
        % Initial position of prey
        disp(' ')
        disp('Select prey centroid')
        [iC.x(2),iC.y(2)] = imInteract(im,'points',1);
        
        % Radius for prey
        disp(' ')
        disp('Select roi radius for prey')
        iC.r(2) = imInteract(im,'radius',iC.x(2),iC.y(2));
        
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
        [iC.xArms,iC.yArms] = imInteract(im,'points',iC.armNum);  
        
        % Radius around mouth
        disp(' ')
        disp('Select radius around body area (without arms)')
        iC.rMouth = imInteract(im,'radius',iC.x,iC.y);
    end
    
    clear im
end


%% Determine area range


if isempty(dir([currDataPath filesep 'Initial conditions.mat']))    
    
    % Frames for plot
    visFrames = round(linspace(iC.startFrame,iC.endFrame,9));
    
    % Make figure
    f = figure;
    
    % Render 9 frames with current settings
    for i = 1:9
        % Get image, subtract around tank
        im{i} = getFrame(vid_path,v,visFrames(i),iC.invert,'gray');
        im{i} = applyMask(im{i},iC.xTank,iC.yTank);
        
        subplot(3,3,i)
        imshow(im{i},'InitialMag','fit')
        title(['Frame ' num2str(visFrames(i))])
    end
    
    brighten(bLevel)
    
    % Index for run thru loop
    runNum = 1;

    % Loop thru for each body
    while true
        
        % Find blob at cX,cY
        [propBod,bwOut] = findBlobs(im{1},imMean,propDiff,'coord advanced',...
                                  iC.x(runNum),iC.y(runNum));
        bodArea = propBod.Area;
                      
        % Default area range                    
        minArea = bodArea/4;
        maxArea = bodArea*10;       

        disp(' ')
        
        if strcmp(iC.expType,'pred prey') && runNum==1
            disp('PREDATOR ----------------------------') 
            clr = pdClr;
            
        elseif strcmp(iC.expType,'pred prey')
            disp('PREY --------------------------------')      
            clr = pyClr;
            
        else
            clr = [0 1 0];
        end
        
        disp('Interactively adjust area bounds');
        disp(' ')
        disp('Up arrow :      Increase max area')
        disp('Down arrow :    Decrease max area')
        disp('Right arrow :   Increase min area')
        disp('Left arrow :    Decrease min area')
        disp(' ')
        
        while true
            
            % Adjust green areas
            for i = 1:9
                % Find blob at cX,cY
                [props,bwOut] = findBlobs(im{i},imMean,propDiff,'all',...
                                          minArea,maxArea);
                
                figure(f);
                
                subplot(3,3,i)
                hold on
                % Make a truecolor all-green image, make non-blobs invisible
                clrField = cat(3, clr(1).*ones(size(im{i})), ...
                                  clr(2).*ones(size(im{i})), ...
                                  clr(3).*ones(size(im{i})));
                h(i) = imshow(clrField,'InitialMag','fit');
                %brighten(bLevel)
                set(h(i), 'AlphaData', bwOut)
                hold off
            end
            
            % Interactive mode
            [x,y,b] = ginput(1);
            
            % Up arrow
            if b==30
                maxArea = round(maxArea + 0.2*maxArea);
                
            % Down arrow
            elseif b==31
                maxArea = max([1 round(maxArea - 0.2*maxArea)]);
                
            % Left arrow
            elseif b==28
                minArea = max([1 round(minArea - 0.2*bodArea)]);
                
            % Right arrow
            elseif b==29
                minArea = round(minArea + 0.2*bodArea);
                
            elseif isempty(b)
                iC.minArea(runNum,1)  = minArea;
                iC.maxArea(runNum,1)  = maxArea;
                break
            end
            
            disp(['Body area = ' num2str(bodArea) ' min = ' num2str(minArea) '  max = ' num2str(maxArea)])
            
            delete(h)
        end
        
        % If first of 2 iterations
        if strcmp(iC.expType,'pred prey') && runNum==1
            runNum = 2;
            
        % Otherwise, leave
        else
            close(f)
            clear frames im
            break
        end
    end
    
    clear runNum minArea maxArea h i x y b clr props bwOut visFrames
end  
   

%% Calibration
    
if isempty(dir([currDataPath filesep 'Initial conditions.mat']))  

    % First image
    im = getFrame(vid_path,v,iC.startFrame,iC.invert,'rgb');
    
    f = figure;
    imshow(im,'InitialMag','fit')
    
     ButtonName = questdlg('Are there calibration landmarks in view?', ...
        '','Yes', 'No', 'Cancel', 'Yes');
    
    if strcmp(ButtonName,'Yes')
        
        % Calibration
        disp(' ')
        disp('Select two points for spatial calibration')
        [xCal,yCal] = imInteract(im,'points',2);
        
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
    
    clear ButtonName xCal yCal answer
    close(f)
    
end

  
    
%% Load/save initial conditions, if present    
    
if ~isempty(dir([currDataPath filesep 'Initial conditions.mat']))  
    load([currDataPath filesep 'Initial conditions'])
else
    % Save initial conditions
    save([currDataPath filesep 'Initial conditions'],'iC')
    
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
    save([currDataPath filesep 'Body'],'Body')
    
elseif strcmp(method,'body arm') && ...
        isempty(dir([currDataPath filesep 'Body arm.mat']))
      
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
    save([currDataPath filesep 'Body arm'],'Body')
    
    % Two body tracking (i.e. predator-prey) ---
elseif strcmp(method,'pred prey') && ...
       isempty(dir([currDataPath filesep 'Body.mat']))
    
    % Make empty Centroid structure
    Body.frames       = [iC.startFrame:iC.endFrame]';
    Body.tform{1,1}   = affine2d(eye(3));
    Body.tform{1,2}   = affine2d(eye(3));
    Body.x            = nan(length(Body.frames),2);
    Body.y            = nan(length(Body.frames),2);
    Body.ang          = nan(length(Body.frames),2);
    %Body.y_flip     = nan(size(Body.frames));
    
    % Initial values
    Body.x(1,:)   = iC.x;
    Body.y(1,:)   = iC.y;
    Body.ang      = [0 0];
    
    % Save data
    save([currDataPath filesep 'Body'],'Body')

end





function clipInfo = selectDuration(vid_path,v)
% Interactively prompts to select a duration for analysis

% Make figure
f = figure;

% Default first and end frames
firstFrame = v.UserData.FirstFrame;
lastFrame  = v.UserData.LastFrame;

% Loop for multiple tries at frame numbers
while true

    % Get images
    %im1 = getFrame(vid_path,v,firstFrame,0,'green screen',imMean);
    %im1 = getFrame(vid_path,v,firstFrame,0,'rgb',imMean);
    im1 = getFrame(vid_path,v,firstFrame,0,'rgb');
    im2 = getFrame(vid_path,v,lastFrame,0,'rgb');
    %im2 = getFrame(vid_path,v,lastFrame,0,'gray',imMean);
    
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



