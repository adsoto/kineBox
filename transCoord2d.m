function varargout = transCoord2d(trans_type,varargin)
% Transforms 2d data from one coordinate system to another

% If data_in are coordinates (n x 2):
%  coord_out =  transCoord2d(trans_type,coord_in,tform);
%  trans_type - type of transformation ('global to local' or 
%               'local to global' for coordinates or 'bw G2L' or 'bw L2G' 
%               for binary images)
%  tform - transformation structure for L system defined in G system, 
%          created by defineSystem2d
% If data_in is an image (n x m):
%  im_out =  transCoord2d(im_,trans_type,tform,im);
%  trans_type - type of transformation ('global to local' or 
%               'local to global' for coordinates or 'bw G2L' or 'bw L2G' 
%               for binary images)
%  tform - transformation structure for L system defined in G system, 
%          created by defineSystem2d
%
% Code developed by McHenryLab at UC Irvine

%% Translate inputs

% Inputs for coordinate transformation
if strcmp(trans_type,'G2L') || ...
   strcmp(trans_type,'L2G')
    
    % First input needs to be tform
    tform = varargin{1};
    
    % Second input needs to be the coordinates in L frame
    coord_in = varargin{2};
 
    % Check coordinate dimensions
    if size(coord_in,2)~=2
        error('Coordinates need to given as a n x 2 matrix')
    end
    
    
% Inputs for image transformation
elseif strcmp(trans_type,'bw L2G')
   
   % First input needs to be tform
    tform = varargin{1};
    
    % Image in roi FOR
    bw_roi = varargin{2};
    
    % binary denoting the roi in global FOR
    bw_G = varargin{3};
    
    % Image in global FOR
    bw_roi_mask = varargin{4};
    
else
    error('trans_type not recognized');
end    

% Check dimensions of tform
if tform.Dimensionality~=2
    error('Code only handles 2D transformations')
end

%% Transformations
    
% Global to local transformation (coordinates)
if strcmp(trans_type,'global to local')
    
    % Translate
    coord_in(:,1) = coord_in(:,1) - tform.T(3,1);
    coord_in(:,2) = coord_in(:,2) - tform.T(3,2);
    
    % Rotate
    coord_out = [tform.T(1:2,1:2) * coord_in']';
    
    
% Local to global transformation (coordinates)
elseif strcmp(trans_type,'local to global')
    
    % Rotate points
    coord_out = (tform.T(1:2,1:2) \ coord_in')';
    
    % Translate global coordinates wrt origin
    coord_out(:,1) = coord_out(:,1) + tform.T(3,1);
    coord_out(:,2) = coord_out(:,2) + tform.T(3,2);
    
% Local to global transformation (coordinates)
elseif strcmp(trans_type,'bw L2G')
    
    % Coordinate system for im_roi
    R = imref2d(size(bw_roi)); 
    
    % Adjust WorldLimits to restrict transformation to just rotation
    % around center
    R.XWorldLimits = R.XWorldLimits-mean(R.XWorldLimits);
    R.YWorldLimits = R.YWorldLimits-mean(R.YWorldLimits);
    
    % Stablize image
    roiRot = imwarp(bw_roi,R,invert(tform),'OutputView',R,...
         'FillValues',255,'SmoothEdges',true);
%     
%     % White out beyond roi
    roiRot(~bw_roi_mask) = 0;
       
    % Black out region outside of roi
    imOut = logical(bw_G.*0);
    
    imOut(bw_G) = roiRot;
    
    
%     % Extract origin
%     origin = tform.T(3,1:2);
%     
%     % Remove translation component to transformation
%     tform.T(3,:) = [0 0 1];
%     
%     % Perform rotation
%     im_rot = imwarp(im_roi,tform,'OutputView',imref2d(size(im_roi)),...
%                   'FillValues',255,'SmoothEdges',true);
%     
%     % Perform translation
%     data_out = imtranslate(im_rot,-origin);
%     
%     
%     aaa=3;
    
% Global to local transformation (coordinates)
elseif strcmp(trans_type,'im G2L')
    
    
else
    
    error('Do not recognize requested transformation')
    
end

%% Set outputs

if strcmp(trans_type,'G2L') || ...
   strcmp(trans_type,'L2G')

    varargout{1} = coord_out;
    
elseif strcmp(trans_type,'bw L2G')
    
    varargout{1}  = imOut;
    
end


% 
% 
% % FUNCTIONS --------------------------
% 
% function ptsT = globalToLocal(tform,coord_in)
% % Assumes column vectors for coordinates
% 
% 
% 
% %pts = [x y];
% 
% % Translate
% coord_in(:,1) = coord_in(:,1) - tform.T(3,1);
% coord_in(:,2) = coord_in(:,2) - tform.T(3,2);
% 
% % Rotate points
% data_out = [tform.T(1:2,1:2) * coord_in']';
% 
% % Extract columns of points
% %xT = ptsT(:,1);
% %yT = ptsT(:,2);
% end
% 
% function ptsT = localToGlobal(tform,pts)
% % Assumes columns vectors for coordinates
% 
% % Check dimensions
% if tform.Dimensionality~=2
%     error('Code only handles 2D transformations')
% end
% 
% % Loop thru columns of coordinates
% i = 1;
%     
%     %pts = [x(:,i) y(:,i)];
%     
%     % Rotate points
%     data_out = (tform.T(1:2,1:2) \ coord_in')';
%     
%     % Translate global coordinates wrt origin
%     data_out(:,1) = data_out(:,1) + tform.T(3,1);
%     data_out(:,2) = data_out(:,2) + tform.T(3,2);
% 
%     
%     clear ptsT pts
% %end
% end