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

if strcmp(bMode,'area')
    
    % Survey blobs
    props = regionprops(bw,'Centroid','Area',...
        'MajorAxisLength','MinorAxisLength',...
        'PixelIdxList','PixelList');
    
    areas = [];
    
    % Loop thru blobs
    for i = 1:length(props)
  
        % If blob is in the area bounds . . .
        if (props(i).Area >= areaMin) && ...
                (props(i).Area <= areaMax)
            
            % Add to area
            areas(j,1) = props(i).Area;
            
            % Add to props out
            propOut(j,1) = props(i);
            
            % Add white pixels for current blob
            bwOut(props(i).PixelIdxList) = 1;
            
            j = j + 1;
            %x = 1;
            %plot(props(i).PixelList(:,1),props(i).PixelList(:,2),'g.')
        end        
    end
    
elseif strcmp(bMode,'coord')
    
    bw = bwselect(bw,x,y);
    
    % Survey blobs
    propOut = regionprops(bw,'Centroid','Area',...
        'MajorAxisLength','MinorAxisLength',...
        'PixelIdxList','PixelList');
    
    if length(props)~=1
        error('Need to select one (and only one) blob')
    end
    
    % Add to area
    areas= props.Area;
            
    % Add white pixels for current blob
    bwOut(props.PixelIdxList) = 1;
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