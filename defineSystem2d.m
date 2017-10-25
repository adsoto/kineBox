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
% S = defineSystem2d('roi',roi0,Centroid,Rotation)
% Defines a region-of-interest within an image
%    roi0 - structure defining region-of-interest for initial frame
%    Centroid - strcuture with centroid data
%    Rotation - (optional) strcuture for rotation data
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
    
    if nargin > 3
        Rotation = varargin{3};
    else
        Rotation = [];
    end
    
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
    
    % Store centroid data
    S.frames        = Centroid.frames;
    S.xCntr         = Centroid.x;
    S.yCntr         = Centroid.y;
    
    % For the 'advanced rotation' method by tracker . . .
    if length(Rotation)==1 && isfield(Rotation,'tform_roi')
        
        S.tform       = Rotation.tform_roi;
        S.ref_frame   = Rotation.ref_frame;
        %S.rot_deg     = Rotation_net;
        
    elseif length(Rotation)~=length(Centroid.x)
        error('mismatch in length of data sources');
        
    end
    
    % Set reference as first frame
    S.ref_frame    = zeros(length(Centroid.x),1);
    S.ref_frame(1) = 1;
    
    % Reference rotation
    ref_rot = 0;
    
    % Loop thru frames, store varying parameters
    for i = 1:length(Centroid.x)
        
        if ~isfield(S,'tform') && ~isempty(Rotation) && length(Rotation)~=1
            S.tform(:,:,i) = Rotation(i).tform_roi;
            
        elseif ~isfield(S,'tform')
            S.tform = [];
        end
        
        % Get angular rotation since last reference frame
        tformInv  = invert(S.tform(:,:,i));
        rot_ang   = atan2(tformInv.T(2,1),tformInv.T(1,1))*180/pi;
            
        %TODO: calculate total rotation
        
        % If this is a reference frame, update reference angle
        if Rotation.ref_frame(i)
            
            % Add current angle to prior reference
            ref_rot = ref_rot + rot_ang;
            
            % Zero rotation angle for storage
            rot_ang = 0;
        end
        
        % Number of roi points
        numroipts = length(roi0.xPerimG);
        
        % Cooridnates of centroid
        xC = Centroid.x(i);
        yC = Centroid.y(i);
        
        % Current roi
        S.roi(i) = giveROI('define','circular',numroipts,roi0.r,xC,yC);
        
        % Store net angle
        S.ang(i,1) = ref_rot + rot_ang;
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


