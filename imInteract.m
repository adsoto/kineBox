function varargout = imInteract(im,action,varargin)
% Image interactive mode 
%   im       - image 
%   action   - strong requesting type of iteraction 
%
% tVal = imInteract(im,'threshold',bLevel)
%    returns the chosen threshold
%    bLevel - brighteness/darkness level (negative values from -1 to 0
%    darken)
%
% [x,y] = imInteract(im,'points',n)
%    returns coordinates for n number of points (default:n=inf)
% r = imInteract(im,'radius')
%    returns radius for a region of interest
% [areaMin,areaMax] = imInteract(im,'area')
%    returns bounds of blob area for undefined threshold
% [areaMin,areaMax] = imInteract(im,'area',tVal)
%    returns bounds of blob area for defined for tVal treshold value
%
% [tVal,x,y] = imInteract(im,'threshold and selection')
%    returns the chosen threshold and blob position
%
% Developed by McHenryLab, UC Irvine


%% Parse inputs

if strcmp(action,'points')
    
    if nargin>2
        n = varargin{1};
    else
        n = inf;
    end
    
    if nargin > 3
        bLevel = varargin{2};
    else
        bLevel = 0;
    end
    
elseif strcmp(action,'radius')
    
    if nargin > 2
        rX = varargin{1};
        rY = varargin{2};
    else
        rX = [];
        rY = [];
    end
    
    if nargin > 4
        bLevel = varargin{3};
    else
        bLevel = 0;
    end
    
elseif strcmp(action,'ellipse')
    
    % Max num of pts
    n = 4;
    
    % Points in each quadrant of ellipse
    numpts = 100;
    
    rX = [];
    rY = [];
    
    if nargin > 2
        bLevel = varargin{1};
    else
        bLevel = 0;
    end
    
elseif strcmp(action,'threshold')
    
    if nargin > 2
        bLevel = varargin{1};
    else
        bLevel = 0;
    end

end


%% Define actions 

% Default interactive mode
iMode = 1;

% Index for the B structure
i = 0;

% Point mode
if strcmp(action,'points')
    % Left click
    i = i + 1;
    B{i}.key = 1;
    B{i}.dostr = 'idx=min([n length(xPos)+1]);xPos(idx,1)=x; yPos(idx,1)=y;';
    B{i}.info = 'Left click: select point';  

    % Right click
    i = i + 1;
    B{i}.key = 3;
    B{i}.dostr = ['if length(xPos)==1;xPos=[];yPos=[]; '...
                  'else; xPos = xPos(1:end-1);yPos = yPos(1:end-1);end'];
    B{i}.info = 'Right click: delete last point';  

% Radius mode
elseif strcmp(action,'radius')
    % Up arrow
    i = i + 1;
    B{i}.key = 30;
    B{i}.dostr = 'r = r + rInc;';
    B{i}.info = 'Up arrow: increase radius';
    
    % Down arrow
    i = i + 1;
    B{i}.key = 31;
    B{i}.dostr = 'r = max([2 (r - rInc)]);';
    B{i}.info = 'Down arrow: decrease radius';    
    
% Ellipse mode
elseif strcmp(action,'ellipse')
    
    disp('Select in clockwise order, starting at the top')
    
    % Left click
    i = i + 1;
    B{i}.key = 1;
    B{i}.dostr = 'idx=min([n length(xPos)+1]);xPos(idx,1)=x; yPos(idx,1)=y;';
    B{i}.info = 'Left click: select point';  

    % Right click
    i = i + 1;
    B{i}.key = 3;
    B{i}.dostr = ['if length(xPos)==1;xPos=[];yPos=[]; '...
                  'else; xPos = xPos(1:end-1);yPos = yPos(1:end-1);end'];
    B{i}.info = 'Right click: delete last point';        

% Threshold mode
elseif strcmp(action,'threshold')
 
    % Up arrow
    i = i + 1;
    B{i}.key = 30;
    B{i}.dostr = 'tVal = min([tVal+0.02 1]);';
    B{i}.info = 'Up arrow: increase threshold';
    
    % Down arrow
    i = i + 1;
    B{i}.key = 31;
    B{i}.dostr = 'tVal = max([tVal-0.02 0]);';
    B{i}.info = 'Down arrow: decrease threshold';
    
% Threshold and selection mode
elseif strcmp(action,'threshold and selection')
 
    % Up arrow
    i = i + 1;
    B{i}.key = 30;
    B{i}.dostr = 'tVal = min([tVal+0.02 1]);';
    B{i}.info = 'Up arrow: increase threshold';
    
    % Down arrow
    i = i + 1;
    B{i}.key = 31;
    B{i}.dostr = 'tVal = max([tVal-0.02 0]);';
    B{i}.info = 'Down arrow: decrease threshold';
    
    % Left click
    i = i + 1;
    B{i}.key = 1;
    B{i}.dostr = 'xPos = x;yPos = y;';
    B{i}.info = 'Left click: select blob';
    
% Area mode
elseif strcmp(action,'area')
    
    % Up arrow
    i = i + 1;
    B{i}.key = 30;
    B{i}.dostr = 'areaMin = areaMin + areaMin_inc;';
    B{i}.info = 'Up arrow: increase min area';
    
    % Down arrow
    i = i + 1;
    B{i}.key = 31;
    B{i}.dostr = 'areaMin = max([0 (areaMin - areaMin_inc)]);';
    B{i}.info = 'Down arrow: decrease min area';
    
    % Right arrow
    i = i + 1;
    B{i}.key = 29;
    B{i}.dostr = 'areaMax = areaMax + areaMax_inc;';
    B{i}.info = '->: increase max area';
    
    % Left arrow
    i = i + 1;
    B{i}.key = 28;
    B{i}.dostr = 'areaMax = max([0 areaMax-areaMax_inc]);';
    B{i}.info = '<-: decrease max area';
    
% Display blobs
elseif strcmp(action,'display blobs')
    % No commands
    B = [];
    
    % Override interactive mode
    iMode = 0;
    
% If no match    
else
    error('Do not recognize action')
end


%% Prep for interaction 

% Figure window
f = figure;

% Give instructions
giveInfo(B)

% Plot image
imshow(im,'InitialMagnification','fit');
brighten(bLevel)
hold on

% Parameter defaults
totArea = size(im,1)*size(im,2);
areaMin = totArea/10^4;
areaMax = totArea/10^2;
areaMin_inc = areaMin/10;
areaMax_inc = areaMax/10;
areaMean = nan;
tVal = graythresh(im);
xPos = [];
yPos = [];
r = round(min(size(im))/4);
rInc = round(min(size(im))/50);

if strcmp(action,'radius') || strcmp(action,'ellipse')
    
    % Radial positions
    theta = linspace(0,2*pi,200);
    
    % Default position
    if isempty(rX)
        rX = size(im,2)/2;
        rY = size(im,1)/2;
    end
 
end

% Overwrite threshold, if provided in Area mode
if strcmp(action,'area')
    % Define threshold, if provided
    if nargin > 2
        tVal = varargin{1};
    end
end

% Overwrite area, if provided in threshold mode
if strcmp(action,'threshold') || ...
   strcmp(action,'threshold and selection')
    areaMin = 0;
    areaMax = inf;
end

% Overwrite area, if provided in threshold mode
if strcmp(action,'display blobs')
    
    % Define threshold and area bounds, if provided
    if nargin > 2
        tVal = varargin{1};       
        if nargin > 3
            areaMin = varargin{2};          
            if nargin > 4
                areaMax =  varargin{3};
            end
        end
    end
end

    
%% Interative mode
    
% Loop interaction
while true
    
    % Show blobs, if needed
    if strcmp(action,'threshold') || strcmp(action,'area')
        % Overlay blobs
        [props,bw,areas,xB,yB] = findBlobs(im,tVal,'area',areaMin,areaMax);
        
        % Make a truecolor all-green image, make non-blobs invisible
        green = cat(3, zeros(size(im)), ones(size(im)), zeros(size(im)));
        h = imshow(green,'InitialMag','fit');
        %brighten(bLevel)
        set(h, 'AlphaData', bw)

    % Show points, if needed
    elseif strcmp(action,'points')
        
        h = plot(xPos,yPos,'g+');
        
    % Show circle, if needed
    elseif strcmp(action,'radius')
        
        % Define circle
        xCirc = r.*cos(theta)+rX;
        yCirc = r.*sin(theta)+rY;
        
        % Plot
        %h = plot(xCirc,yCirc,'g-');
        h = line(xCirc,yCirc,'Color',[0 1 0 0.2],'LineWidth',3);
        
    % Ellipse mode
    elseif strcmp(action,'ellipse')
        
        % Radial positions
        theta1 = linspace(-pi/2,0,numpts)';
        theta2 = linspace(0,pi/2,numpts)';
        theta3 = linspace(pi/2,pi,numpts)';
        theta4 = linspace(pi,1.5*pi,numpts)';
        
        if length(xPos)<2
            h = plot(xPos,yPos,'g+');
            
        else
            radY1 = abs(diff(yPos(1:2)));
            radX1 = abs(diff(xPos(1:2)));
            
            % Define circle
            xCirc = radX1.*cos(theta1)+xPos(1);
            yCirc = radY1.*sin(theta1)+yPos(2);
            
            if length(xPos)>2
                
                radY2 = abs(diff(yPos(2:3)));
                radX2 = abs(diff(xPos(2:3)));
                
                % Define circle
                xCirc = [xCirc; radX2.*cos(theta2)+xPos(3)];
                yCirc = [yCirc; radY2.*sin(theta2)+yPos(2)];
                
                if length(xPos)>3
                    
                    radY3 = abs(diff(yPos(3:4)));
                    radX3 = abs(diff(xPos(3:4)));
                    
                    % Define circle
                    xCirc = [xCirc; radX3.*cos(theta3)+xPos(3)];
                    yCirc = [yCirc; radY3.*sin(theta3)+yPos(4)];
                    
                    radY4 = abs(yPos(4)-yPos(1));
                    radX4 = abs(xPos(4)-xPos(1));
                    
                    % Define circle
                    xCirc = [xCirc; radX4.*cos(theta4)+xPos(1)];
                    yCirc = [yCirc; radY4.*sin(theta4)+yPos(4)];
                end
            end
                
            h = plot(xCirc,yCirc,'g-');
            
        end
            
%         elseif length(xPos)==3
%             
%             % Radial positions
%             theta = linspace(-pi/2,pi/2,200);
%             
%             radY =  range(yPos)/2;
%             radX =  range(xPos);
%             
%             % Define circle
%             xCirc = radX.*cos(theta) + mean([xPos(1) xPos(3)]);
%             yCirc = radY.*sin(theta) + mean(yPos);
%             
%             % Define circle
%             xCirc = r.*cos(theta)+xPos(1);
%             yCirc = r.*sin(theta)+mean(yPos);
%             
%             h = plot(xCirc,yCirc,'g-');
%             
%         elseif length(xPos)>3
%             
%             % Radial positions
%             theta = linspace(0,2*pi,200);
%             
%             radY =  range(yPos);
%             radX =  range(xPos);
%             
%             % Define circle
%             xCirc = radX.*cos(theta) + mean(xPos(1));
%             yCirc = radY.*sin(theta) + mean(yPos);
%             
%             h = plot(xCirc,yCirc,'g-');    
%         end
        
    end

    % Tile: threshold mode
    if strcmp(action,'threshold')
        title(['threshold = ' num2str(round(tVal*255),'%4.0f')])
        
    % Tile: area mode    
    elseif strcmp(action,'area')
        title(['A_m_i_n = ' num2str(round(areaMin),'%4.0f') ...
            ', A_m_a_x = ' num2str(round(areaMax),'%4.0f') ...
            ', A_m_e_a_n = ' num2str(round(mean(areas)),'%4.0f') ])   
        
    elseif strcmp(action,'points')
        title(['Num points = ' num2str(length(xPos),'%4.0f')])  
        
    elseif strcmp(action,'radius')
        title(['r = ' num2str(r,'%4.2f') ' pix'])      
        
    end
     
    % Interacive mode 1 (response to single input)
    if iMode==1
        % Get input
        [x,y,b] = ginput(1);
        
        % If return pressed
        if isempty(b)
            break
        end
        
        % Loop thru keystroke commands
        for i = 1:length(B)
            if b==B{i}.key
                eval(B{i}.dostr);
                break
            end
        end
        
    % Break loop, if not in interactive mode
    else     
        hold off
        break
    end
    
    % Remove perimeter
    delete(h)
end

% Close figure
close(f);


%% Selection processing

if strcmp(action,'threshold and selection')       
    % Choose blob
    bw = bwselect(bw,xPos,yPos);
    
    % Get properties of blobs
    props = regionprops(bw,'Centroid','Area',...
                        'MajorAxisLength','MinorAxisLength');
                    
    % Check that one selected               
    if length(props)>1    
        error('More than one object selected');
        
    elseif isempty(props)
        error('No object selected');
    end   
end


%% Define output

% Threshold mode
if strcmp(action,'threshold')
    varargout{1} = tVal;
    
elseif strcmp(action,'area')
    varargout{1} = areaMin;
    varargout{2} = areaMax;
    
elseif strcmp(action,'threshold and selection')
    varargout{1} = tVal;
    varargout{2} = xPos;
    varargout{3} = yPos;
    
elseif strcmp(action,'points')
    varargout{1} = xPos;
    varargout{2} = yPos;
    
elseif strcmp(action,'radius')
    varargout{1} = r;
    
elseif strcmp(action,'ellipse')
    varargout{1} = xCirc;
    varargout{2} = yCirc;
end



function giveInfo(B)
% Report instructions for each keystroke
if ~isempty(B)
    disp(' ')
    disp('Press keys for the following actions.')
    for i = 1:length(B)
        disp(['    ' B{i}.info])
    end
    disp('Press return when finished')
    disp(' ')
end



