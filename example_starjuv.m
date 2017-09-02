function example_starjuv

%% Parameters

% Factor by which the blob diameter gets multiplied for roi
blobMultiple = 1.3;


%% Preliminaries

% Extract root directories
paths = givePaths;

% Sample for analysis
vid_path = [paths.vid_root filesep 'Seastars' filesep 'CSULB juvenile' ...
    filesep 'SS18' filesep 'timelapse7'];
        
% Load video info (v)
v = defineVidObject(vid_path,'JPG');

[x,y] = tracker(vid_path,v,'threshold',1.3,1);

% return
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






