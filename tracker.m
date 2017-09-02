function varargout = tracker(vid_path,v,method,roiScale,visTracking,frames)
% Tracks a landmark or body
%
%  vid_path - path to video file or image sequence
%  v - structure of info about video (generated by defineVidObject)
%  method - type of tracking approach: 'threshold' (currently only method)
%  roiScale - Factor by which the blob diameter gets multiplied for roi
%  visTracking - Logical to visualize the tracking
%  frames - listing of frame numbers to analyze
%
% [x,y] = threshTrack - returns just centroid cooridnates
%
% Developed by McHenryLab at UC Irvine


%% Parameter defaults

% If frames undefined
if nargin < 6
    % Frame numbers
    frames = v.UserData.FirstFrame:v.UserData.LastFrame;
    
    % If visTracking undefined . . .
    if nargin < 5
        % Turn on
        visTracking = 1;
        
        % If roiScale undefined . . .
        if nargin<4
            % Factor by which the blob diameter gets multiplied for roi
            roiScale = 1.3;
            
            % If method undefined . . .
            if nargin<3
                method = 'threshold';
            end
        end
    end
end

% Check method input
if ~strcmp(method,'threshold')
    error('requested method not recognized');
end

% Angular coordinates for circular roi
theta = linspace(0,2*pi,500);


%% Prompt for user input on first frame

% First image
im = rgb2gray(getFrame(vid_path,v,frames(1)));

f = figure('DoubleBuffer','on');

% Select sea star from first frame
[tVal,x,y] = imInteract(im,'threshold and selection');

% If no visualization . . .
if ~visTracking
    % Pause and close figure
    pause(0.5)
    close
else
    hold off
end


%% Find blob diameter

% Define binary image of seastar
bw = im2bw(im,tVal);
bw = bwselect(bw,x,y);
bw = imfill(bw,'holes');

% Find perimeter
[yPerim, xPerim] = find(bwperim(bw,8));

% Blob diameter
blobDiam = max([range(xPerim) range(yPerim)]);

clear xPerim yPerim bw


%% Tracking object

% Loop thru frames
for i = 1:length(frames)
    
    % Current frame
    cFrame = frames(i);
    
    % Current image
    im = rgb2gray(getFrame(vid_path,v,cFrame));

    % Make blob image 
    bw = im2bw(im,tVal);
    bw = bwselect(bw,x,y);
    bw = imfill(bw,'holes');
    
    % Find center of blob
    props = regionprops(bw,'Centroid');
    
    % Check for one (and only one) blob        
    if length(props)>1    
        error('More than one object in roi');       
    elseif isempty(props)
        error('No object in roi');
    end
    
    % Circular coordinates for new roi 
    r     = roiScale*blobDiam/2;
    xC    = r.*cos(theta) + props.Centroid(1);
    yC    = r.*sin(theta) + props.Centroid(2);

    % Visualize
    if visTracking
        imshow(im,'InitialMagnification','fit');
        hold on
        h = line(xC,yC,'Color',[0 1 0 0.2],'LineWidth',3);
        title(['Frame ' num2str(cFrame) '/' num2str(frames(end))])
        plot(props.Centroid(1),props.Centroid(2),'g+')
        hold off
        pause(0.01)    
    else
        disp(['Done tracking frame ' num2str(cFrame) '/' num2str(frames(end))])
    end
    
    % Store centroid
    x(i,1) = props.Centroid(1);
    y(i,1) = props.Centroid(2);
    
    % Clear for next loop
    clear cFrame r xC yC props bw bwMask im cFrame
end


%% Define outputs

% If fewer than 2 outputs
if nargout<2
    error('There should at least be two outputs to this function')

% Otherwise
else
    varargout{1} = x;
    varargout{2} = y;
    
    % If more than 2
    if nargout>2
        varargout{3} = f;
    end
end

if ~visTracking
    disp(' ')
    disp('Trackign complete!')
end
