function varargout = findBlobs(im,imMean,propDiff,bMode,varargin)
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
    
    if length(varargin) > 2
        specialAction = varargin{3};
    else
        specialAction = [];
    end
    
    if length(varargin) > 3
        minArea = varargin{4};
    else
        minArea = 0;
    end
    
    if length(varargin) > 4
        maxArea = varargin{5};
    else
        maxArea = inf;
    end   
    
elseif strcmp(bMode,'advanced comparison')    

     % Initial blob
    prop0 = varargin{1};

    % Previous position
    x = varargin{2};
    y = varargin{3};
    
    
    % Special action
    specialAction = varargin{4};
    
    % Area bounds
    areaMin = varargin{5};
    areaMax = varargin{6};
    
%      % areaFactor - proportional range in area wrt last blob
%     if length(varargin) > 3
%         
%     else
%         areaFactor = 0.5;
%     end
    
    
%     x = propLast.Centroid(1);
%     y = propLast.Centroid(2);
    
    
elseif strcmp(bMode,'area and circ')
    
    areaMin = varargin{1};
    areaMax = varargin{2};
    AR_max  = varargin{3};
    
elseif strcmp(bMode,'coord')
    x = varargin{1};
    y = varargin{2};
    
elseif strcmp(bMode,'all')
    
    if length(varargin)>0
        areaMin = varargin{1};
    else
        areaMin = 0;
    end
    
    if length(varargin)>1
        areaMax = varargin{2};
    else
        areaMax = inf;
    end
       
else
    error('bMode not recognized');
end


%% Make binary image

% If color image . . .
if size(im,3)==3
    % Convert images to HSV
    H      = rgb2hsv(im);
    Hmean  = rgb2hsv(imMean);
    
    % Make binary, based on differences from mean image in hue and value
    bw = ( H(:,:,1)>(1+propDiff*2)*Hmean(:,:,1)) ...
        | (H(:,:,1)< (1-propDiff*2) * Hmean(:,:,1)) ...
        & H(:,:,3)<(0.99)*Hmean(:,:,3);
    
% If grayscale . . .
else
    bw = im < (1-propDiff)*imMean;
end

clear H Hmean

% Fill holes
bw = imfill(bw,'holes');

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
    
elseif strcmp(bMode,'all')
    
    % Survey blobs
    props = regionprops(bw,'Centroid','Area',...
        'MajorAxisLength','MinorAxisLength',...
        'PixelIdxList','PixelList');
      
    j = 1;
    for i = 1:length(props)
        if (props(i).Area > areaMin) && ...
                (props(i).Area < areaMax)
            
            propOut(j,1) = props(i);
            j = j + 1;
        end
    end
    
    if j == 1
        propOut = [];
    else
        propOut = props;
        
        % Add to area
        areas = propOut.Area;
        
        for j = 1:length(propOut)
            % Add white pixels for current blob
            bwOut(propOut(j).PixelIdxList) = 1;
        end
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
    
    % Trim fins, if requested
    if strcmp(specialAction,'trim fins')
        bw2 = bwdist(~bw);
        bw2 = imadjust(bw2./max(bw2(:)));
        tVal = min([1 graythresh(bw2)]);
        bw = im2bw(bw2,tVal);
        
        clear bw2 tVal
    end
    
    % Dialate & erode the binary image a bit
    se = strel('disk',3,4);    
    bw = imdilate(bw,se);
    bw = imerode(bw,se);
    
    % Try to get blob on coordinate ------------
    bwS = bwselect(bw,x,y);
    
    % Survey blobs
    propOut = regionprops(bwS,'Centroid','Area',...
        'MajorAxisLength','MinorAxisLength',...
        'PixelIdxList','PixelList');
    
    % If that fails, decide among blobs
    if isempty(propOut)
        
        % Survey blobs
        props = regionprops(bw,'Centroid','Area',...
            'MajorAxisLength','MinorAxisLength',...
            'PixelIdxList','PixelList');
        
        minDist = inf;
        
        if isempty(props)
            error('Lost blob -- maybe try different treshold or roi radius')
          
        elseif length(props)==1
            propOut = props;
            
        elseif length(props)>1 
            for i = 1:length(props)
                
                % Distance of current blob from last
                currDist = hypot(x-props(i).Centroid(1),y-props(i).Centroid(2));
                
                if (props(i).Area > minArea) && ...
                        (currDist < minDist)  && ...
                        (props(i).Area < maxArea)
                    
                    minDist = currDist;
                    propOut = props(i);
                end
            end
            
            % If lost blob, try again without area filter
            if isempty(propOut)
                for i = 1:length(props)
                    
                    % Distance of current blob from last
                    currDist = hypot(x-props(i).Centroid(1),y-props(i).Centroid(2));
                    
                    if currDist < minDist
                        minDist = currDist;
                        propOut = props(i);
                    end
                end 
            end
            
        end
    end

    % Add to area
    areas = propOut.Area;
            
    % Add white pixels for current blob
    bwOut(propOut.PixelIdxList) = 1;

    
elseif strcmp(bMode,'advanced comparison')
    
    % Trim fins, if requested
    if strcmp(specialAction,'trim fins')
        bw2 = bwdist(~bw);
        bw = bw2>max(bw2(:)*.25);
%         bw2 = imadjust(bw2./max(bw2(:)));
%         tVal = min([1 2*graythresh(bw2)]);
%         bw = im2bw(bw2,tVal);
        
        clear bw2 tVal
    end
    
    % Dialate & erode the binary image a bit
    se = strel('disk',3,4);    
    bw = imdilate(bw,se);
    bw = imerode(bw,se);
  
    % Survey blobs
    props = regionprops(bw,'Centroid','Area',...
        'MajorAxisLength','MinorAxisLength',...
        'PixelIdxList','PixelList');

    if isempty(props)
        warning('Lost blob');
        propOut = nan;
        
    elseif length(props)==1
        propOut = props;
        
    elseif length(props)>1 
        
        % Initialize index
        k = 1;
        
        % Loop thru blobs
        for j = 1:length(props)
            
            % Include only those within area range
            if props(j).Area > areaMin && ...
               props(j).Area < areaMax
           
                % Distance from last
                dCenter(k) = hypot(x-props(j).Centroid(1),y-props(j).Centroid(2));
                
                % Store blob
                props2(k) = props(j);
                
                k = k + 1;
            end
        end
        
        % If blobs passed area filter, choose closest from last
        if k > 1
            
            % Best is closest
            iBest = find(dCenter==min(dCenter),1,'first');
            
            % Advance best
            propOut = props2(iBest);
        
        % If no blobs passed filter, advance nan
        else
            propOut = nan;
            
        end
    end

    if isstruct(propOut)
        % Add white pixels for current blob
        bwOut(propOut.PixelIdxList) = 1;
    end
    
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