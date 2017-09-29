function varargout = tracker(vid_path,v,method,roiScale,visTracking,frames, ...
                             Centroid,Rotation,makeVid)
% Tracks a landmark or body
%
%  vid_path - path to video file or image sequence
%  v - structure of info about video (generated by defineVidObject)
%  method - type of tracking approach: 'threshold' (currently only method)
%  roiScale - Factor by which the blob diameter gets multiplied for roi
%  visTracking - Logical to visualize the tracking
%  frames - listing of frame numbers to analyze
%  Centroid - structure with fields x_pix, y_pix, & r_pix for roi coord & radius
%
% [x,y] = threshTrack - returns just centroid cooridnates
%
% Developed by McHenryLab at UC Irvine


%% Parameter defaults

% If frames undefined, or empty
if nargin < 6 || isempty(frames)
    
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

% If threshold method . . .
if strcmp(method,'threshold')

    % Define frames
    Centroid.frames = frames;   
    
% If body rotation . . .
elseif strcmp(method,'body rotation')
    
    % Make sure Centroid data provided
    if nargin < 7 
        error(['Centroid structure needs to be defined before analyzing '...
               'body rotation. Try threshold method to resolve centroid values']);
    end
   
% If visualization . . .
elseif  strcmp(method,'visualize') 
    
    % Make figure
    fig = figure;
    
    % If no rotation data
    if nargin < 8
        
        includeRot = 0;
        
        if nargin < 7
            error(['Centroid structure needs to be defined before analyzing '...
                'body rotation. Try threshold method to resolve centroid values']);
        end
        
    % If there is rotation data
    else
        includeRot = 1;
    end
    
% If no match on method
else
    error('requested method not recognized');
end

% Angular coordinates for circular roi
theta = linspace(0,2*pi,500);


%% Prompt for user input on first frame

% For threshold method
if strcmp(method,'threshold')
    % First image
    im = rgb2gray(getFrame(vid_path,v,frames(1)));
    
    f = figure('DoubleBuffer','on');
    
    % Select sea star from first frame
    [tVal,x,y] = imInteract(im,'threshold and selection');
    
    Centroid.x_pix(1) = x;
    Centroid.y_pix(1) = y;
    
    % If no visualization . . .
    if ~visTracking
        % Pause and close figure
        pause(0.5)
        close
    else
        hold off
    end
    
    clear x y 
end


%% Initial parameter values

% For threshold method: find initial blob diameter & perimeter
if strcmp(method,'threshold')
    
    % Define binary image of seastar
    bw = im2bw(im,tVal);
    bw = bwselect(bw,Centroid.x_pix,Centroid.y_pix);
    bw = imfill(bw,'holes');
    
    % Find perimeter
    [yPerim, xPerim] = find(bwperim(bw,8));
    
    % Blob diameter & radius
    blobDiam = max([range(xPerim) range(yPerim)]);
    r        = roiScale*blobDiam/2;
    
    % Store radius
    Centroid.r_pix = r;
    
    clear xPerim yPerim bw r
    
elseif strcmp(method,'body rotation')
    
    % Initialize image registration parameters
    [optimizer, metric]  = imregconfig('monomodal');
    optimizer.MaximumStepLength = 5e-4;
    optimizer.MaximumIterations = 1500;
    optimizer.RelaxationFactor  = 0.2;
    
    % First image
    im0 = rgb2gray(getFrame(vid_path,v,frames(1)));
    
    % Focus on roi
    im0 = isolate_roi(im0,Centroid.x_pix(1),Centroid.y_pix(1),Centroid.r_pix,theta);
end


%% Tracking object

% Loop thru frames
for i = 1:length(frames)
    
    % Current frame
    cFrame = frames(i);
    
    % Current image
    im = rgb2gray(getFrame(vid_path,v,cFrame));
    
    % Threshold method: find centroid coordinates
    if strcmp(method,'threshold')
        
        % Find centroid for current frame
        [tmp_x,tmp_y] = findCentroid(im,tVal,Centroid.r_pix,theta,visTracking,...
            Centroid.x_pix(end),Centroid.y_pix(end),cFrame,frames);
        
        % Store results
        Centroid.x_pix(i,1) = tmp_x;
        Centroid.y_pix(i,1) = tmp_y;
        
        % Don't include rotation in plotting
        includeRot = 0;
        
        % Clear for next loop
        clear tmp_x tmp_y
        
    % body rotation method: 
    elseif strcmp(method,'body rotation')
        
        tform = findRot(im,im0,Centroid.x_pix(i),Centroid.y_pix(i),...
                        Centroid.r_pix,theta,optimizer,metric);
          
        %disp(['Rotation: completed ' num2str(i) ' of ' num2str(length(frames))])            
         % Include rotation in plotting
         includeRot = 1;
%         
%         % Focus on roi
%         [im_roi,imMask] = isolate_roi(im,Centroid.x_pix(i),Centroid.y_pix(i), ...
%                                       Centroid.r_pix,theta);
%         
%         % Transformation object to stablize head wrt imHead0
%         tform_roi = imregtform(im_roi,im0,'rigid',optimizer,metric);
%          
%         % Stablize head image
%         imStable = imwarp(im_roi,tform_roi,'OutputView',imref2d(size(im0)));
%         
%         % White out beyond roi
%         imStable(~imMask) = 255;
        
        
        
%         if visTracking
%             
%             % Circular coordinates for new roi
%             xC    = Centroid.r_pix.*cos(theta) + Centroid.x_pix(i);
%             yC    = Centroid.r_pix.*sin(theta) + Centroid.y_pix(i);
%           
%             subplot(1,2,1)
%             imshow(im,'InitialMagnification','fit');
%             hold on
%             h = line(xC,yC,'Color',[0 1 0 0.2],'LineWidth',3);
%             title(['Frame ' num2str(cFrame) '/' num2str(frames(end))])
%             plot(Centroid.x_pix(i),Centroid.y_pix(i),'g+')
%             hold off
%             
%             subplot(1,2,2)
%             imshow(imStable,'InitialMagnification','fit');
%             
%             pause(0.01)
            
%         else
            
            % Status
            disp(['Rotation: completed ' num2str(i) ' of ' num2str(length(frames))])
            
%         end

        % Store results
        Rotation(i).tform_roi = tform;
        
        % Clear for next iteration
        clear im_roi tform_roi imStable xC yC h
    end
    
    % Visualize: 
    if strcmp(method,'visualize') || visTracking
        
        % Title text
        t_txt = ['Frame ' num2str(cFrame) '/' num2str(frames(end))];
            
        % If rotation data included . . .
        if includeRot
            visTrack(im,Centroid.x_pix(i),Centroid.y_pix(i),Centroid.r_pix,theta,...
                Rotation(i).tform_roi,t_txt);
            
        % If no rotation
        else

            % Visualization code
            visTrack(im,Centroid.x_pix(i),Centroid.y_pix(i),Centroid.r_pix,...
                theta,[],t_txt);
        end
        
        % Log frame
        if makeVid
            M(i) = getframe(gcf);
        end
        
        % Pause briefly to render
        pause(0.001)
    end

end


%% Define outputs

% Threshold method
if strcmp(method,'threshold')
    
        varargout{1} = Centroid;
        
    
% Body rotation method    
elseif strcmp(method,'body rotation')
    
    varargout{1} = Rotation;
    

elseif strcmp(method,'visualize') 
    
    if makeVid
        varargout{1} = M;
    else
        varargout{1} = fig;
    end
end

% Report status
if ~visTracking
    disp(' ')
    disp([method ' complete!'])
end


function tform_roi = findRot(im,im0,x,y,r,theta,optimizer,metric)
% Find rotation matrix using image registration

% Focus on roi
[im_roi,imMask] = isolate_roi(im,x,y,r,theta);

% Transformation object to stablize head wrt im0
tform_roi = imregtform(im_roi,im0,'rigid',optimizer,metric);

% if visTracking
%     
%     % Stablize image
%     %imStable = imwarp(im_roi,tform_roi,'OutputView',imref2d(size(im0)));
%     imStable = imwarp(im_roi,tform_roi);
%     
%     % White out beyond roi
%     imStable(~imMask) = 255;
%     
%     % Circular coordinates for new roi
%     xC    = Centroid.r_pix.*cos(theta) + Centroid.x_pix(i);
%     yC    = Centroid.r_pix.*sin(theta) + Centroid.y_pix(i);
%     
%     subplot(1,2,1)
%     imshow(im,'InitialMagnification','fit');
%     hold on
%     h = line(xC,yC,'Color',[0 1 0 0.2],'LineWidth',3);
%     title(['Frame ' num2str(cFrame) '/' num2str(frames(end))])
%     plot(Centroid.x_pix(i),Centroid.y_pix(i),'g+')
%     hold off
%     
%     subplot(1,2,2)
%     imshow(imStable,'InitialMagnification','fit');
%     
%     pause(0.01)
%     
% else
    
    % Status
    
    
%end

function visTrack(im,x,y,r,theta,tform,t_txt)

% If rotation included . . 
if nargin>5 && ~isempty(tform)
    % Focus on roi
    [im_roi,imMask] = isolate_roi(im,x,y,r,theta);
    
    % Stablize image
    imStable = imwarp(im_roi,tform,'OutputView',imref2d(size(im_roi)));
    
    % White out beyond roi
    imStable(~imMask) = 255;

    subplot(1,2,2)
    imshow(imStable,'InitialMagnification','fit');
    
    % Set up for main plot
    subplot(1,2,1)
end

% Circular coordinates for new roi
xC    = r.*cos(theta) + x;
yC    = r.*sin(theta) + y;

imshow(im,'InitialMagnification','fit');
hold on
h = line(xC,yC,'Color',[0 1 0 0.2],'LineWidth',3);
title(t_txt)
plot(x,y,'g+')
hold off


function [x,y] = findCentroid(im,tVal,r,theta,visTracking,xLast,yLast,cFrame,frames)
% Used by threshold method to locate the center of a blob

% Make blob image
bw = im2bw(im,tVal);
bw = bwselect(bw,xLast,yLast);
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
xC    = r.*cos(theta) + props.Centroid(1);
yC    = r.*sin(theta) + props.Centroid(2);

% Visualize
if visTracking
%     imshow(im,'InitialMagnification','fit');
%     hold on
%     h = line(xC,yC,'Color',[0 1 0 0.2],'LineWidth',3);
%     title(['Frame ' num2str(cFrame) '/' num2str(frames(end))])
%     plot(props.Centroid(1),props.Centroid(2),'g+')
%     hold off
%     pause(0.01)
else
    disp(['Done centroid frame ' num2str(cFrame) '/' num2str(frames(end))])
end

% Store centroid
x = props.Centroid(1);
y = props.Centroid(2);


function [im,roi_mask] = isolate_roi(im,x,y,r,theta)

% Maximum size of an image dimension
maxSize = 250;

% % rectangular ROI vector
% roi_rect = [Centroid.x_pix(1)-Centroid.r_pix ...
%     Centroid.y_pix(1)-Centroid.r_pix ...
%     Centroid.r_pix*2 ...
%     Centroid.r_pix*2];

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


