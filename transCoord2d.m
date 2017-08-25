function coord_out = transCoord2d(coord_in,tform,trans_type)
% Transforms 2d coordinates from one system to another
%  coord_in - input coordinates (n x 2)
%  tform - transformation structure for L system defined in G system, 
%          created by defineSystem2d
%  trans_type - type of tranformation ('global to local' or 
%               'local to global')
%
% Code developed by McHenryLab at UC Irvine


% Check coordinate dimensions
if size(coord_in,2)~=2
    error('Coordinates need to given as a n x 2 matrix')
end

% Check dimensions of tform
if tform.Dimensionality~=2
    error('Code only handles 2D transformations')
end


% Global to local transformation
if strcmp(trans_type,'global to local')
    
    % Translate
    coord_in(:,1) = coord_in(:,1) - tform.T(3,1);
    coord_in(:,2) = coord_in(:,2) - tform.T(3,2);
    
    % Rotate
    coord_out = [tform.T(1:2,1:2) * coord_in']';
    
    
% Local to global transformation
elseif strcmp(trans_type,'local to global')
    
    % Rotate points
    coord_out = (tform.T(1:2,1:2) \ coord_in')';
    
    % Translate global coordinates wrt origin
    coord_out(:,1) = coord_out(:,1) + tform.T(3,1);
    coord_out(:,2) = coord_out(:,2) + tform.T(3,2);
    
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
% coord_out = [tform.T(1:2,1:2) * coord_in']';
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
%     coord_out = (tform.T(1:2,1:2) \ coord_in')';
%     
%     % Translate global coordinates wrt origin
%     coord_out(:,1) = coord_out(:,1) + tform.T(3,1);
%     coord_out(:,2) = coord_out(:,2) + tform.T(3,2);
% 
%     
%     clear ptsT pts
% %end
% end