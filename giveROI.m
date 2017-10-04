
function varargout = giveROI(roi_type,varargin)
% Returns data and images for a region-of-interest (ROI)
%
%[bw_mask,im_roi,roi_rect,bw_roi_mask] = giveROI(...
%   bw_mask - binary image for the roi with dimensions of im
%   im_roi - image within the roi (downsampled, if dSampling==1)
%   roi_rect - [x y width height] of the roi bounding box
%   bw_roi_mask - binary of mask in roi
%
% ... = giveROI('circular',im,xC,yC,r,theta,dSampling)
%   im - image
%   xC - roi x-center
%   yC - roi y-center
%   r - roi circle radius
%   theta - angular position
%   dSampling - logical that indicates whether to downsample roi image
%
% Developed by McHenryLab at UC Irvine


%% Translate inputs

if strcmp(roi_type,'circular')
    
    im         = varargin{1};
    xC         = varargin{2};
    yC         = varargin{3};
    r          = varargin{4};
    theta      = varargin{5};
    dSampling  = varargin{6};
    
else
    error('Do not recognize roi_type')
end


%% Parameters

% Maximum size of an image dimension (for downsampling)
maxSize = 250;

% rectangular ROI vector
roi_rect = [round(xC-r) round(yC-r) ceil(r*2) ceil(r*2)];

% Circular coordinates for new roi in image FOR
xCirc    = round(floor(r).*cos(theta)+ xC);
yCirc    = round(floor(r).*sin(theta) + yC);

% Circular coordinates for new roi in roi FOR
xCirc_roi    = round(floor(r).*cos(theta) + r);
yCirc_roi    = round(floor(r).*sin(theta) + r);


%% Create images

% Binary mask in image im
bw_mask = roipoly(size(im,1),size(im,2),xCirc,yCirc);

% Crop image
im_roi  = imcrop(im,roi_rect);

% White out pixels outside of circular roi
bw_roi_mask = roipoly(size(im_roi,1),size(im_roi,2),yCirc_roi,xCirc_roi);

% White out area around roi
im_roi(~bw_roi_mask) = 255;

% Reduce size (helps speed up image registration)
if dSampling && (length(im)>maxSize)
    
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


%% Define output

varargout{1} = bw_mask;
varargout{2} = im_roi;
varargout{3} = roi_rect;
varargout{4} = bw_roi_mask;


