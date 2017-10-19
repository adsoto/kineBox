function varargout = findBlobs(im,tVal,bMode,varargin)
% Finds blobs defined by threshold and either area bounds or coordinates
%
% bMode = 'area' : findBlobs(im,tVal,bMode,areaMin,areaMax)
% bMode = 'coord' : findBlobs(im,tVal,bMode,x,y)
% [props,bwOut] = findBlobs(im,tVal,...) - returns props structure for blobs                                           
% [props,bwOut,areas] = findBlobs(im,tVal,...) - also returns blob areas
% [props,bwOut,areas,x,y] = findBlobs(im,tVal,...) - also returns
%                                                   peripheral coordinates
%


%% Handle inputs

% Set mode-specific parameters
if strcmp(bMode,'area')
    areaMin = varargin{1};
    areaMax = varargin{2};

elseif strcmp(bMode,'coord advanced')
    x        = varargin{1};
    y        = varargin{2};
    if nargin<5
        areaLast = varargin{3};
    else
        areaLast = [];
    end
    
elseif strcmp(bMode,'area and circ')
    areaMin = varargin{1};
    areaMax = varargin{2};
    AR_max  = varargin{3};
    
elseif strcmp(bMode,'coord')
    x = varargin{1};
    y = varargin{2};
    
else
    error('bMode not recognized');
end


%% Make binary image

% Make binary
bw = im2bw(im,tVal);

% Fill holes
bw = imfill(~bw,'holes');

% Start with black image
bwOut = bw.*0~=0;   


%% Select blobs, according to mode

% Initialize index
j = 1;

if strcmp(bMode,'area') || strcmp(bMode,'area and circ')
    
    % Survey blobs
    props = regionprops(bw,'Centroid','Area',...
        'MajorAxisLength','MinorAxisLength',...
        'PixelIdxList','PixelList');
    
    areas = []; 
    
    % Loop thru blobs
    for i = 1:length(props)
  
        % If blob is in the area bounds . . .
        if (props(i).Area >= areaMin) && (props(i).Area <= areaMax)
            
            if strcmp(bMode,'area')
                % Add to area
                areas(j,1) = props(i).Area;
                
                % Add to props out
                propOut(j,1) = props(i);
                
                % Add white pixels for current blob
                bwOut(props(i).PixelIdxList) = 1;
                
                j = j + 1;
            
            elseif strcmp(bMode,'area and circ') && ...
                    props(i).MajorAxisLength/props(i).MinorAxisLength < AR_max
                
                % Add to props out
                propOut(j,1) = props(i);
                
                % Add white pixels for current blob
                bwOut(props(i).PixelIdxList) = 1;
                
                j = j + 1;
            end
            
            %x = 1;
            %plot(props(i).PixelList(:,1),props(i).PixelList(:,2),'g.')
        end        
    end
    
    % If no blobs, return empty variable
    if j==1
        propOut = [];
    end
    
elseif strcmp(bMode,'coord')
    
    bw = bwselect(bw,x,y);
    
    % Survey blobs
    propOut = regionprops(bw,'Centroid','Area',...
        'MajorAxisLength','MinorAxisLength',...
        'PixelIdxList','PixelList');
    
    if length(propOut)~=1
        error('Need to select one (and only one) blob')
    end
    
    % Add to area
    areas = propOut.Area;
            
    % Add white pixels for current blob
    bwOut(propOut.PixelIdxList) = 1;
    
    
elseif strcmp(bMode,'coord advanced')
    
    % Dialate the binary image a bit
    se = strel('disk',3,4);    
    bw = imdilate(bw,se);
    
    % Try to get blob on coordinate ------------
    bwS = bwselect(bw,x,y);
    
    % Survey blobs
    propOut = regionprops(bwS,'Centroid','Area',...
        'MajorAxisLength','MinorAxisLength',...
        'PixelIdxList','PixelList');
    
    % If that fails, get blob closest to last
    if isempty(propOut)
        
        % Survey blobs
        props = regionprops(bw,'Centroid','Area',...
            'MajorAxisLength','MinorAxisLength',...
            'PixelIdxList','PixelList');
        
        minDist = [inf nan];
        
        if isempty(props)
            error('Lost blob -- maybe try different treshold or roi radius')
            
        elseif length(props)>1
            for j = 1:length(props)  
                
                % Distance of current blob from last
                currDist = hypot(x-props(i).Centroid(1),y-props(i).Centroid(2));
                
                if ~isempty(areaLast) && ...
                   (props(i).Area > areaLast*.5) && ...
                   (currDist < minDist)     
                    minDist = currDist;
                    propOut = props(i);
                end
            end
  
        end
    end   
    
    if isempty(propOut)
        error('Lost blob -- maybe try different treshold or roi radius')
    end
    
    % Add to area
    areas = propOut.Area;
            
    % Add white pixels for current blob
    bwOut(propOut.PixelIdxList) = 1;

end

 % Trace perimeter
[yOut, xOut] = find(bwperim(bwOut,8));
        

%% Outputs

% Define outputs
varargout{1} = propOut;
varargout{2} = bwOut;

% Areas, if requested
if nargout>2
    
    varargout{3} = areas;

    if nargout==4
        error('This function is organized to return 3 or 5 outputs, not 4');
    end
    
    if nargout > 4
        varargout{4} = xOut;
        varargout{5} = yOut;
    end
    
end

% Overlay
%h = plot(x,y,'.g');

%meanArea = mean(areas);