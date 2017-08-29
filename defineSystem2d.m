function tform = defineSystem2d(origin,axCoord,axType)
% Defines coordinate local system (L) within global system (G)
%    origin  - a vector of 2 coordinates that define origin L in G
%    axCoord - a vector of 2 coordinates defining an axis for L in G
%    axType  - defines which axis of L given by axVect ('x' or 'y')
%
% Code developed by McHenryLab at UC Irvine


%% Check inputs

% Check dimensions of origin
if length(origin)~=2 
    error('Origin needs to have a length of 2')
end

% Check dimensions of axCoord
if length(axCoord)~=2
    error('axCoord needs to have a length of 2')
end

% Check input of axType for 2D
if ~strcmp(axType,'x') && ~strcmp(axType,'y') && ~strcmp(axType,'z')
    error('axType needs to be given as either "x", "y", or "z"')
end


%% Adjust dimensions

% Make origin a column vector
if size(origin,1)>size(origin,2)
    origin = origin';
end

% Make axCoord a column vector
if size(axCoord,1)>size(axCoord,2)
    axCoord = axCoord';
end

%% Translate wrt origin

axCoord(1) = axCoord(1) - origin(1);
axCoord(2) = axCoord(2) - origin(2);


%% Define system from x-axis coordinate

if strcmp(axType,'x')
    
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

if strcmp(axType,'y')
    
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


%% Package system

% Create rotation matrix (from inertial axes to local axes)
R = [xaxis; yaxis; [origin 1]];

% Format trans matrix for matlab
tform = affine2d(R);

