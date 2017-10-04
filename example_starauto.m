function example_starauto
% Sample code for running automatic kinematics of sea star tube feet
% Working with Andres' blue led images


%% Preliminaries

% Extract directories
paths = givePaths;

frInterval = 500;

% Number of frames used for mean image
numFrames = 50;

numThresh = 4;

% % Path to sample zebrafish video
% vid_path = [paths.vid_root filesep 'Seastars' filesep 'Star Prints' ...
%             filesep 'C0423'];

% Path to sample zebrafish video
vid_path = [paths.vid_root filesep 'Seastars' filesep 'Star Prints' ...
            filesep 'C0423'];


% Load video info (v)
v = defineVidObject(vid_path);


%% Prompt user for input

% Read next image
im = getFrame(vid_path,v,1);

% Convert to grayscale, enhance contrast
%imGray = imadjust(adapthisteq(rgb2gray(im)));

% Fine approximate threshold value
%tVal = graythresh(imGray);
%thresh = multithresh(imGray,4);
%imSeg = imquantize(imGray,thresh);
%RGB = label2rgb(imSeg);

%imshow(RGB);
% Binary image
%imBW = im2bw(imGray,tVal);


frSkip = floor(frInterval/numFrames) - 1;

% Loop thru frames
for i = ceil(frInterval/2):(v.UserData.NumFrames-floor(frInterval/2))
    
    frames = (ceil(i-frInterval/2)+1):frSkip:floor(i+frInterval/2);
    
    im = rgb2gray(getFrame(vid_path,v,i));
    
    imMean = meanImage(vid_path,v,'enhance contrast',frames);
    
    %im1 = imsubtract(im,imMean);
    %imMean = adapthisteq(stasisImage(vid_path,v,'enhance contrast',frames));
    
    thresh = multithresh(imMean,numThresh);
    
    imSeg = imquantize(imMean,thresh);
    %RGB = label2rgb(imSeg);
    
    
    
    for j = 1:numThresh
        imFeet = showFeet(imSeg==j,[30 120]);
        
    end
    
%    [centers, radii] = imfindcircles(adapthisteq(imMean),[50 120],'ObjectPolarity','bright', ...
%    'Sensitivity',.9);
    
    

    if 1
        imshow(imMean,'InitialMagnification','fit');
        plot(centers(:,1),centers(:,2),'+r')
        %h = viscircles(centers,radii);
        sss=2
    end
    
    % Update status
    disp(['Done ' num2str(i) ' of ' ...
          num2str((v.UserData.NumFrames-frInterval)) ' frames']);
end

% Display
figure
subplot(1,2,1)
imshow(im,'InitialMagnification','fit');
subplot(1,2,2)
imshow(imGray,'InitialMagnification','fit');

delete(v)




function imOut = showFeet(im,dia_range)

%propDiff

minArea = 0.2*min(dia_range)^2;
minArea2 = min(dia_range)^2;
maxArea = max(dia_range)^2;

se = strel('disk',round(min(dia_range)/4));


im = bwareaopen(im,minArea);

im1 = imdilate(im,se);

im2 = imfill(im1,'holes');

%se2 = strel('disk',round(min(dia_range)/8));

im3 = imerode(im2,se);

im4 = bwareaopen(im3,minArea2);



% Get area and centroid
    props = regionprops(logical(im4),'Centroid','Area',...
                        'MajorAxisLength','MinorAxisLength');
    
    cCent = [];

    for i = 1:length(props)        
        if props(i).Area < maxArea && ...
                props(i).MajorAxisLength
            % Circle centroid
            cCent = [cCent;props(i).Centroid(1) props(i).Centroid(2)];
        end
    end
    
    % Choose largest blob
    if ~isempty(cCent)
        im5 = bwselect(im4,cCent(:,1),cCent(:,2));
    else
        im5 = im4;
    end

imOut = im5;



