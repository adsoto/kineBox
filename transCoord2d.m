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
elseif strcmp(trans_type,'im G2L') || ...
       strcmp(trans_type,'im L2G')
   
   % First input needs to be tform
    tform = varargin{1};
    
    % First input needs to be tform
    im_roi = varargin{2};
    
    % First input needs to be tform
    im_G = varargin{3};
    
else
    error('trans_type not recognized');
end    

% Check dimensions of tform
if tform.Dimensionality~=2
    error('Code only handles 2D transformations')
end
    
% Global to local transformation (coordinates)
if strcmp(trans_type,'global to local')
    
    % Translate
    coord_in(:,1) = coord_in(:,1) - tform.T(3,1);
    coord_in(:,2) = coord_in(:,2) - tform.T(3,2);
    
    % Rotate
    data_out = [tform.T(1:2,1:2) * coord_in']';
    
    
% Local to global transformation (coordinates)
elseif strcmp(trans_type,'local to global')
    
    % Rotate points
    data_out = (tform.T(1:2,1:2) \ coord_in')';
    
    % Translate global coordinates wrt origin
    data_out(:,1) = data_out(:,1) + tform.T(3,1);
    data_out(:,2) = data_out(:,2) + tform.T(3,2);
    
% Local to global transformation (coordinates)
elseif strcmp(trans_type,'im L2G')
    
    % Extract origin
    origin = tform.T(3,1:2);
    
    % Remove translation component to transformation
    tform.T(3,:) = [0 0 1];
    
    % Perform rotation
    im_rot = imwarp(im_roi,tform,'OutputView',imref2d(size(im_roi)),...
                  'FillValues',255,'SmoothEdges',true);
    
    % Perform translation
    data_out = imtranslate(im_rot,-origin);
    
    
    aaa=3;
    
% Global to local transformation (coordinates)
elseif strcmp(trans_type,'im G2L')
    
    
else
    
    error('Do not recognize requested transformation')
    
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