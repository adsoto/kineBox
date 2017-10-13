
function varargout = giveROI(action,varargin)
% Returns features for a region-of-interest (ROI)
% giveROI(action)
%  action - particular operation relating to a ROI ('define')
%
% roi = giveROI('define',shape)
% Defines a region of interest
%   roi - structure of roi properties 
%   shape - string indicating the shape of the roi ('circular')
%
% roi = giveROI('define','circular',num_pts,r)
%   num_pts - number of points for bounding the roi
%   r - radius of circular roi
%   Centroid - structure definining the center coordinates for the roi
%
% roi = giveROI('define','circular',num_pts,r,xC,yC)
%   xC - roi x-center
%   yC - roi y-center
%
% im_roi = giveROI('unstabilized',im,roi,dSample)
% returns roi image, without rotation correction
%   im - full frame image
%   roi - structure with stats on roi
%   dSample - logical that indicates whether to downsample the roi image
%
% im_roi = giveROI('stabilized',im,roi,dSample)
% returns roi image, WITH rotation correction
%
%[im_roi,bw_mask,roi_rect,bw_roi_mask,imStable] = giveROI(...
%   im_roi - image within the roi (downsampled, if dSampling==1)
%   bw_mask - binary image for the roi with dimensions of im
%   roi_rect - [x y width height] of the roi bounding box
%   bw_roi_mask - binary of mask in roi
%   imStable - image stablized against rotation
%
% Developed by McHenryLab at UC Irvine


%% Translate inputs 

if strcmp(action,'define')
    
    roi.shape       = varargin{1};
    num_pts         = varargin{2};
    roi.r           = varargin{3};
    
    % Define center points
    if nargin>4
        roi.xCntr = varargin{4};
        roi.yCntr = varargin{5};
    else
        roi.xCntr = 0;
        roi.yCntr = 0;
    end
    
    
elseif strcmp(action,'unstabilized')
    
    im      = varargin{1};
    roi     = varargin{2};
    dSample = varargin{3};
%     xC         = varargin{2};
%     yC         = varargin{3};
%     r          = varargin{4};
%     theta      = varargin{5};
%     
%     if nargin >6
%         dSampling  = varargin{6};
%     else
%         dSampling = 0;
%     end
%     
%     if nargin >7
%         tform = varargin{7};
%         
%         if nargout<5
%             error(['You need to assign at least 5 outputs for imStable ' ...
%                  'when you provide a tform'])
%         end
%     else
%         tform.T = [];
%     end
    
elseif strcmp(action,'stabilized')
    
    im      = varargin{1};
    roi     = varargin{2};
    dSample = varargin{3};
    tform   = varargin{4};
    
else
    error('Do not recognize action')
end


%% Parameters

if strcmp(action,'define')
    
    % Angular coordinates for circular roi
    roi.theta = linspace(0,2*pi,num_pts);
    
    % Circular coordinates for new roi in image FOR
    roi.xPerimG  = round(floor(roi.r-1).*cos(roi.theta) + roi.xCntr);
    roi.yPerimG  = round(floor(roi.r-1).*sin(roi.theta) + roi.yCntr);
    
    % Circular coordinates for new roi in roi FOR
    roi.xPerimL    = round(floor(roi.r-1).*cos(roi.theta) + roi.r);
    roi.yPerimL    = round(floor(roi.r-1).*sin(roi.theta) + roi.r);
    
    % Bounding rectange for the roi
    roi.rect = [round(roi.xCntr-roi.r) round(roi.yCntr-roi.r) ceil(roi.r*2) ceil(roi.r*2)];

else
    
    % Maximum size of an image dimension (for downsampling)
    maxSize = 250;
%     
%     % rectangular ROI vector
%     if nargin < 9
%         roi_rect = [round(xC-r) round(yC-r) ceil(r*2) ceil(r*2)];
%     else
%         roi_rect = varargin{8};
%     end
%     
%     % Circular coordinates for new roi in image FOR
%     xCirc    = round(floor(r).*cos(theta)+ xC);
%     yCirc    = round(floor(r).*sin(theta) + yC);
%     
%     % Circular coordinates for new roi in roi FOR
%     xCirc_roi    = round(floor(r).*cos(theta) + r);
%     yCirc_roi    = round(floor(r).*sin(theta) + r);
%     
end

% %% Defining an roi
% 
% if strcmp(action,'define')
%     
%     
% end


%% Create images

if ~strcmp(action,'define')
    % Binary mask in image im
    %bw_mask = roipoly(size(im,1),size(im,2),xCirc,yCirc);
%     bw_mask = roipoly(size(im,1),size(im,2),...
%         [roi_rect(1) roi_rect(1)+roi_rect(3)+1 ...
%         roi_rect(1)+roi_rect(3)+1 roi_rect(1) ...
%         roi_rect(1)],...
%         [roi_rect(2) roi_rect(2) ...
%         roi_rect(2)+roi_rect(4)+1 roi_rect(2)+roi_rect(4)+1 ...
%         roi_rect(2)]);
    
    bw_mask = roipoly(size(im,1),size(im,2),...
        [roi.rect(1) roi.rect(1)+roi.rect(3)+1 ...
         roi.rect(1)+roi.rect(3)+1 roi.rect(1) ...
         roi.rect(1)],...
        [roi.rect(2) roi.rect(2) ...
        roi.rect(2)+roi.rect(4)+1 roi.rect(2)+roi.rect(4)+1 ...
        roi.rect(2)]);
    
    % Crop image
    im_roi  = imcrop(im,roi.rect);
    
    % White out pixels outside of circular roi
    %bw_roi_mask = roipoly(size(im_roi,1),size(im_roi,2),yCirc_roi,xCirc_roi);
    bw_roi_mask = roipoly(size(im_roi,1),size(im_roi,2),roi.yPerimL,roi.xPerimL);
    
    % White out area around roi
    im_roi(~bw_roi_mask) = 255;
    
    % Reduce size (helps speed up image registration)
    if dSample && (length(im)>maxSize)
        
        % Factor by which to resize
        imFactor = maxSize/length(im);
        
        % Downsample image
        im_roi = imresize(im_roi,imFactor);
        
        bw_roi_mask = imresize(bw_roi_mask,imFactor);
        % White out pixels outside of circular roi
        %bw_roi = roipoly(size(im_roi,1),size(im_roi,2),round(yC*imFactor),round(xC*imFactor));
    end
    
    % Check that roi is square
    if size(im_roi,1)~=size(im_roi,2)
        error('Image not square and needs to be');
    end
    
    if strcmp(action,'stabilized')
        
        if isempty(tform.T)
            error('You need to provide tform if you want imStable');
        end
        
        % Coordinate system for im_roi
        R = imref2d(size(im_roi));
        
        % Adjust WorldLimits to restrict transformation to just rotation
        % around center
        R.XWorldLimits = R.XWorldLimits-mean(R.XWorldLimits);
        R.YWorldLimits = R.YWorldLimits-mean(R.YWorldLimits);
        
        % Stablize image
        imStable = imwarp(im_roi,R,tform,'OutputView',R,...
            'FillValues',255,'SmoothEdges',true);
        
        % White out beyond roi
        imStable(~bw_roi_mask) = 255;
        
        % Deliver imStable
        im_roi = imStable;
    end
end


%% Define output

% If defining roi . . .
if strcmp(action,'define')
    
    varargout{1} = roi;
    
% Otherwise . . .
else
    
    varargout{1} = im_roi;
    
    if nargout>1
        varargout{2} = bw_mask;
        
        if nargout>2
            varargout{3} = roi_rect;
            
            if nargout>3
                varargout{4} = bw_roi_mask;
                
            end
        end
    end
    
    
    
end


