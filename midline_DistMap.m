function blob = midline_DistMap(blob)
% MIDLINE_DISTMAP Finds the midline of a fish using distance mapping
% INPUTS:   
%           blob        = image blob
% 
% Note that this function operates on a single frame, the function should
% be called within a loop that iterates through a sequence of frames. 
%
% Developed by McHenryLab at UC Irvine

%% Inputs and parameters

% Indicator to visualize progress
visProg = 1;

% If no input blob given, generate one
if nargin < 1
    
    % indicator to use static background as mean image
    newMean = 1;
    
    % Prompt to browse to filename
    [filename,pathname] = uigetfile({'*.jpg';'*.tif'});

%     % If there is a video info dir . . .
%     if ~isempty(dir([pathname filesep 'vid_info.mat']))
%         
%         % Load video info (v)
%         load([pathname filesep 'vid_info'])
%     else
%         % Warning message, 
%         warning('Video info missing, run makeImageSeq to proceed')
%         
%         % Prompt user to run makeImageSeq
%         prompt = 'Would you like to run makeImageSeq now? Y or N  ';
%         str = input(prompt,'s');
%         
%         if strcmp('Y',str) || strcmp('y',str)
%         
%             % Run makeImageSeq
%             makeImageSeq
%             
%             % Load video info (v)
%             load([pathname filesep 'vid_info'])
%             
%         elseif strcmp('N',str) || strcmp('n',str)
%             return
%         end      
% end
    
    % Read first frame in image sequence folder
    [im,~] = imread([pathname filesep filename]);
        
    % Get background image
    imMean = makeMeanImage(pathname,newMean);
    
    % Find blobs that define the fish in frame
    blob = findBlobs(im,imMean,0);
end

% Size of binary image (j rows, k columns)
szBlob = size(blob.BW);

%% FIND TAIL POINT 

% Get centroid of small blob
stats = regionprops(blob.BWsmall,'Centroid','Area');

% Skeletonize blob by thinning
skel = bwmorph(blob.BW,'thin',Inf);

% Find endpoints
endPts = bwmorph(skel, 'endpoints');

% Get coordinates for end points
[yEnd,xEnd] = find(endPts);

% Distance between centroid and endpoints
xCent = repmat(stats.Centroid(1),length(xEnd),1);
yCent = repmat(stats.Centroid(2),length(yEnd),1);
xDiff = xEnd-xCent;
yDiff = yEnd-yCent;
dists = hypot(xDiff,yDiff);

% Index for new tail point as max distance to endpoint
[~,iVal] = max(dists);
xTail = xEnd(iVal);
yTail = yEnd(iVal);

% Clear extraneous variables
clear endPts xEnd yEnd xCent yCent dists iVal k D

%% Distance map along rows

% Generate distance map with 'euclidean' option
D1 = bwdist(~blob.BW,'euclidean');

% Find maximum pixel distance along each row 
% NOTE: the j_th entry in colsIdx gives column index for max value of row j
[maxRows,colsIdx] = max(D1,[],2); 

% Find maximum of all pixels (by looking through the max of each row)
[maxPix1,colIdx] = max(maxRows);

% Rows to exclude (whenever distance value is zero)
exCol = maxRows<1;

% Vector for plotting row results: 1:#rows
nRows = (1:szBlob(1))';

%% Distance map along columns

% Find maximum pixel distance along each column
% NOTE: the k_th entry in rowsIdx gives row index for max value of column k
[maxCols,rowsIdx] = max(D1,[],1); 

% Find maximum of all pixels (by looking through the max of each column)
[maxPix2,rowIdx] = max(maxCols);

% Columns to exclude (whenever distance value is zero)
exRow = maxCols<1;

% Vector for plotting column results: 1:#columns
nCols = 1:szBlob(2);

% Plots for debugging
if visProg
    % Convert distance map to RGB for plotting
    RGB1 = repmat(mat2gray(D1), [1 1 3]);
    figure, imshow(RGB1)
    
    % Countour lines (level set) of ditance map
    hold on, imcontour(D1)
    
    % Max distance results from looking along rows
    hold on, plot(colsIdx(~exCol),nRows(~exCol),'+')
    
    % Max distance results from looking along columns
    hold on, plot(nCols(~exRow),rowsIdx(~exRow),'+')
    
    % Plot absolute max value (point near head)
    hold on, plot(rowIdx,colIdx,'go')
    nCols;
end

%% Work with merged data

% Create background image, same size as blob.im
midTrace = false(size(blob.im));

% Vertical concatenation of column numbers to keep for spline: x-values
xData = [colsIdx(~exCol);(nCols(~exRow))'];

% Horizontal concatenation of row numbers to keep for spline: y-values
yData = [(nRows(~exCol))',rowsIdx(~exRow)];

% Convert subscript indices to linear indices
lin_indx = sub2ind(size(midTrace), yData', xData);

% Add points from initial skeleton for full midline
lin_indx = [lin_indx; find(skel)];

% Get rid of repeats
lin_indx = unique(lin_indx,'stable');

% Make midline points white (set to 1)
midTrace(lin_indx) = 1;

% Structuring element for morphological operations
se = strel('disk',2);

% Dilate midline points (helps to connect any disjoint regions)
midDilate = imdilate(midTrace,se);

% Skeletonize midline blob & get rid of small blobs
skel = bwmorph(midDilate,'thin',Inf);
skel = bwareaopen(skel, 15);

% Find endpoint coordinates in new skeleton 
[yEnd,xEnd] = find(bwmorph(skel,'endpoints'));

% Compute distance between endpoints and previous tail point
distTail = hypot((xTail-xEnd),(yTail - yEnd));

% Index of boundary point nearest to previous tail point
[~,tInd] = min(distTail);

% Update tail point
xTail = xEnd(tInd);
yTail = yEnd(tInd);

% Compute distance map from tail with 'bwdistgeodesic' & 'chessboard' option
% Note: 'bwdistgeodesic' takes in the mask of the skeleton which specifies
% which pixels the path is allowed to traverse
tailD = bwdistgeodesic(skel,xTail,yTail,'chessboard');

% Set NaNs to 0
tailD(isnan(tailD)) = 0;

% Find maximum pixel distance along each column of tailD
[maxD,rowsIdxD] = max(tailD,[],1);

% Find maximum of all pixels (by looking through the max of each column)
[~,headIdx] = max(maxD);

% Row of point furthest from tail
yHead = rowsIdxD(headIdx);

% Column of point furthest from tail
xHead = headIdx;

% Recompute distance map from tail using 'bwdistgeodesic'
tailD = bwdistgeodesic(skel,xTail,yTail,'quasi-euclidean');

% Compute distance map from head using 'bwdistgeodesic'
headD = bwdistgeodesic(skel,xHead,yHead,'quasi-euclidean');

% Add tail and head distance maps and round the result
totD = tailD + headD;
totD = round(totD * 16) / 16;

% Convert NaNs into Inf
totD(isnan(totD)) = inf;

% Find the shortest path along the skeleton
skeleton_path = imregionalmin(totD);

% Get coordinates of midline
[yPoints,xPoints] = find(skeleton_path);

% Plot shortest path along midline for debugging
if 0
    % pad solution path for plotting
    thick_sol_path = imdilate(skeleton_path, ones(3,3));
    
    % Create figure window
    figure(15), 
    
    % Overlay the pixels on the shortest path in midTrace in yellow.
    P = imoverlay(midTrace, thick_sol_path, [1 1 0]);
    imshow(P, 'InitialMagnification', 'fit')
end

% Convert subscript indices of refined midline to linear indices
lin_indx2= sub2ind(size(headD), yPoints, xPoints);

% Distance values from head of refined midline points
dist_mid = headD(lin_indx2);

% Sort distance values (and keep track of positions)
[sort_dist_mid,sort_indx] = sort(dist_mid);

% Sorted x-coordinate of midline
blob.xMid = xPoints(sort_indx);

% Sorted y-coordinate of midline
blob.yMid = yPoints(sort_indx);

% Arclength along midline
blob.sMid = sort_dist_mid;

% Visualize progress up to this point
if visProg
    % Create figure window
    figure(20)
    
    % Show current image frame
    imshow(blob.im,[],'InitialMagnification','fit')
    hold on
    
    % Overlay midline points
    plot(blob.xMid,blob.yMid,'g-')
    hold off
%     pause(0.1)
end

function blob = findBlobs(imStart,imMean,adjustON)
% FINDBLOBS takes in the current image, performs a background subtraction
% and finds any blobs in the current frame
%
% INPUTS: imStart   = current frame
%         imMean    = background or mean image for background subtraction
%         adjustON  = indicator variable for automatic contrast correction
%

% Minimum blob size in pixels, smaller objects will be discarded
minBlob = 15;

% Adjust grayscale values if adjustON indicator is set to 1
if adjustON
    im     = (imadjust(imStart));
    imSub  = (imadjust(imMean));
else
    im = imStart;
    imSub = imMean;
end

% Subtract background
warning off
im = imsubtract(imSub,im);
warning on

% Get inverse of image
im = imcomplement(im);

% Find threshold
tVal = min([0.95 graythresh(im)+0.1]);

% Threshold image
imBW    = ~im2bw(im,tVal);

% Close image to get rid of small contamination 
se    = strel('diamond',3);
imBW = imclose(imBW,se);

% Get peripheral shapes, query size and boundaries
blobShapes = regionprops(imBW,'ConvexArea','BoundingBox',...
                  'Orientation','Centroid','Image','FilledImage');

% If no fish or touching a wall . . .
if isempty(blobShapes) 
    % nan blob
    imBlob = nan;
    
    % Store nans
    blob.xPerim   = nan;
    blob.yPerim   = nan;
    blob.roi_blob = nan;
    
% If fish . . .
else
    
    % Filter out small objects (< minBlob px) found by 'regionprops'
    blobShapes = blobShapes([blobShapes.ConvexArea] > minBlob);
    
    % Find largest object and get is bounding box
    [~,ind2] = max([blobShapes.ConvexArea]);
    perim2 = blobShapes(ind2).BoundingBox;
    
    % Number of pixels that pad the blob
    pad_val = 10;
    
    % Rectangle for blob roi: [XMIN YMIN WIDTH HEIGHT]
    rect = [perim2(1) - pad_val,    perim2(2) - pad_val, ...
            perim2(3) + 2*pad_val,  perim2(4) + 2*pad_val];
                     
    % Crop down image to predator region using 'BoundingBox' output
    % crop image syntax: imcrop(im,[XMIN YMIN WIDTH HEIGHT])
    imBlob = imcrop(im,rect);
    
    % Crop down binary image to predator region 
    imBlobBW = imcrop(imBW,rect);
    
    % Store predator roi data
    blob.roi_blob  = rect;
    
end

% If there is a blob, find binary blobs (big and small)
if ~isnan(blob.roi_blob(1))
    
    % Find largest blob
    [imBW,~] = returnBlob(imBlobBW);
    
    % Create tmp image (uses blob as mask)
    tmp = imBlob.*0+2.^16;
    tmp(imBW(:)) = imBlob(imBW(:));
    
    % Starting low threshold value (for small blob)
    tVal2 = tVal.*.75;
    
    % Increment to increase threshold (for small blob)
    tIncr = .05;
    
    % Loop to create image at lowest threshold that finds single blob
    while true
        
        % Find blob
        [imBWsmall,numBlob] = returnBlob(~im2bw(tmp,tVal2));
        
        % Check number of blobs
        if numBlob==1
            break
        else
            tVal2 = tVal2 + tIncr;
        end
        
        % Check for too high a threshold
        if tVal2 >= tVal
            tVal2 = tVal;
            imBWsmall = imBW;
            break
        end
    end

% Otherwise, return nans
else
    imBW      = nan;
    imBWsmall = nan;
end

% Store images
blob.im      = imBlob;          % cropped image of fish
blob.BW      = imBW;            % binary image of fish
blob.BWsmall = imBWsmall;       % binary image of fish w/ smaller threshold


function [imBW,numBlob] = returnBlob(imBW)
% Returns binary image of blobs
%
% INPUTS: imBW = cropped binary image of fish blob
%

% Close gaps with dilation and erosion 
se   = strel('disk',5);
imBW = imdilate(imBW,se);
imBW = imerode(imBW,se);

% Identify blobs (BWCONNCOMP is more memory efficient than BWLABEL) 
CC = bwconncomp(imBW);
LL = labelmatrix(CC);

% Get peripheral boundary of shapes
[bb,~] = bwboundaries(imBW,'noholes');

% Select blob with greatest periphery
maxB = 0;
idx = [];

for j = 1:length(bb)
    if length(bb{j}) > maxB
        maxB = length(bb{j});
        idx = j;
    end
end

if isempty(idx)
    imBW = nan;
    numBlob = 0;
else
    % Define image as having only largest blob
    imBW = LL==idx;
    
    % Return number of blobs
    numBlob = max(LL(:));
end