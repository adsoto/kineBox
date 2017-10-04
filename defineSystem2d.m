function tform = defineSystem2d(inType,varargin)
% Defines coordinate local system (L) within global system (G)
%    inType - type of data used to define coordinate system
%    ('x-axis','y-axis')
%
% tform = defineSystem2d('roi tform',rect,tform)
%    rect - vector defining the roi in G frame
%    tform - transformation matrix for rotation within the roi
%
% tform = defineSystem2d('x-axis',origin,axCoord)
%    origin  - a vector of 2 coordinates that define origin L in G
%    axCoord - a vector of 2 coordinates defining x-axis for L in G
%
% tform = defineSystem2d('y-axis',origin,axCoord)
%    origin  - a vector of 2 coordinates that define origin L in G
%    axCoord - a vector of 2 coordinates defining y-axis for L in G
%
% Code developed by McHenryLab at UC Irvine


%% Translate inputs

% If using axis coordinates
if strcmp(inType,'x-axis') || strcmp(inType,'y-axis')
    
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

    
    
elseif strcmp(inType,'roi tform')    
    % Region of interest rectangle
    rect = varargin{1};
    
    % Set axis
    tform = varargin{2};
    
else
    error('inType not recignized');
end


%% Define system from x-axis coordinate

if strcmp(inType,'x-axis')
    
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

if strcmp(inType,'y-axis')
    
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

if strcmp(inType,'roi tform')    
   
    % Redefine origin, leave the rotation matrix
    tform.T(3,:) = [rect(1) rect(2) 1];
    
    % Coordinates for roi
    %tform.roi = rect;
    
end

%% Package system for output

if strcmp(inType,'x-axis') || strcmp(inType,'y-axis')
    % Create rotation matrix (from inertial axes to local axes)
    R = [xaxis; yaxis; [origin 1]];
       
    % Format trans matrix for matlab
    tform = affine2d(R);
end


