function varargout = defineSystem2d(coordType,varargin)
% Defines coordinate local system (L) within global system (G)
% S = defineSystem2d(coordType)
%    S - structure that defines a coordinate system
%    coordType - type of data used to define coordinate system
%    ('x-axis','y-axis','roi')
%
% S = defineSystem2d('x-axis',origin,axCoord)
%    origin  - a vector of 2 coordinates that define origin L in G
%    axCoord - a vector of 2 coordinates defining x-axis for L in G
%
% S = defineSystem2d('y-axis',origin,axCoord)
%    origin  - a vector of 2 coordinates that define origin L in G
%    axCoord - a vector of 2 coordinates defining y-axis for L in G
%
% S = defineSystem2d('full roi',rect,tform)
% Defines a region-of-interest within an image
%    rect - vector defining the roi in G frame
%    tform - transformation matrix for rotation within the roi
%
% Code developed by McHenryLab at UC Irvine


%% Translate inputs

% If using axis coordinates
if strcmp(coordType,'x-axis') || strcmp(coordType,'y-axis')
    
    % Set origin
    origin = varargin{1};
    
    % Set axis
    axCoord = varargin{2};
    
    % Check dimensions of origin
    if length(origin)~=2
        error('Origin needs to have a length of 2')
    end
    
    % Check dimensions of axCoord
    if length(axCoord)~=2
        error('axCoord needs to have a length of 2')
    end
    
    % Make origin a column vector
    if size(origin,1)>size(origin,2)
        origin = origin';
    end
    
    % Make axCoord a column vector
    if size(axCoord,1)>size(axCoord,2)
        axCoord = axCoord';
    end
    
    % Translate wrt origin
    axCoord(1) = axCoord(1) - origin(1);
    axCoord(2) = axCoord(2) - origin(2);
    
elseif strcmp(coordType,'roi')    
    
    % Region of interest rectangle
    roi0 = varargin{1};
    
    % 
    Centroid = varargin{2};
    
    Rotation = varargin{3};
    
else
    error('coordType not recognized');
end


%% Define system from x-axis coordinate

if strcmp(coordType,'x-axis')
    
    % Define xaxis
    xaxis = axCoord;
    
    % Normalize to create a unit vector
    xaxis = [[xaxis./norm(xaxis)] 0];
    
    % Define y-axis
    yaxis = [-xaxis(2) xaxis(1) 0];
    
    % Normalize to create a unit vector
    yaxis = yaxis./norm(yaxis);
     
    % Define z-axis
    zaxis = cross(xaxis,yaxis);
end


%% Define system from y-axis coordinate

if strcmp(coordType,'y-axis')
    
    % Define xaxis
    yaxis = axCoord;
    
    % Normalize to create a unit vector
    yaxis = [[yaxis./norm(yaxis)] 0];
    
    % Define x-axis
    xaxis = [yaxis(2) -yaxis(1) 0];
    
    % Normalize to create a unit vector
    yaxis = yaxis./norm(yaxis);
   
    % Define y-axis
    xaxis = [-xaxis(2) xaxis(1) 0];
 
    % Define z-axis
    zaxis = cross(xaxis,yaxis);
end


%% Define system for an roi

if strcmp(coordType,'roi')
    
    if length(Rotation)~=length(Centroid.x)
        error('mismatch in length of data sources');
    end
    
    % Store general parameters
    S.frames        = Centroid.frames;
    S.xCntr         = Centroid.x;
    S.yCntr         = Centroid.y;
%    S.roi.r         = roi0.r;
%     S.roi_shape     = roi0.shape;
%     S.roi.xPerimL   = roi0.xPerimL;
%     S.roi.yPerimL   = roi0.yPerimL;
%     S.roi.theta     = roi0.theta;
    
    % Loop thru frames, store varying parameters
    for i = 1:length(Rotation)
        
        S.tform(:,:,i) = Rotation(i).tform_roi;
        
        % Number of roi points
        numroipts = length(roi0.xPerimG);
        
        % Cooridnates of centroid
        xC = Centroid.x(i);
        yC = Centroid.y(i);
        
        % Current roi
        S.roi(i) = giveROI('define','circular',numroipts,roi0.r,xC,yC);
%         
%         % Store roi data
%         
%         S.roi(i).xCntr    = Centroid.x(i);
%         S.roi(i).yCntr    = Centroid.y(i);
%         S.roi(i).xPerimG  = tmp.xPerimG;
%         S.roi(i).yPerimG  = tmp.yPerimG;
%         S.roi(i).         = tmp.rect;       
    end
end


%% Package system for output

if strcmp(coordType,'x-axis') || strcmp(coordType,'y-axis')
    % Create rotation matrix (from inertial axes to local axes)
    R = [xaxis; yaxis; [origin 1]];
       
    % Format trans matrix for matlab
    S.tform = affine2d(R);
end

% Output
varargout{1} = S;


