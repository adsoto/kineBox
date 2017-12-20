function [x,y] = videoGUI(vid_path,v,frames,imInvert,acqMode,varargin)
% Interactive acquisition of coordinate points from a video 
% 

%% Parse inputs

if nargin< 5
    acqMode = 'simple';
end

if strcmp(acqMode,'simple')
    
    % roi radius
    if length(varargin)>0
        r = varargin{1};
    else
        r = 0;
    end
    
    % Marker color
    if length(varargin)>1
        mClr = varargin{2};
    else
        mClr = [0 1 0];
    end
    
    % Starting values
    if length(varargin)>3
        xStart = varargin{3};
        yStart = varargin{4};
    else
        xStart = nan(length(frames),1);
        yStart = nan(length(frames),1);
    end
else
    error('Do not recognize acqMode');
end


%% Default parameters

% Container of handle data
H = [];

% Adjust window docking
set(0,'DefaultFigureWindowStyle','normal')

 % Initial frame
im0 = getFrame(vid_path,v,frames(1),imInvert,'gray');    

% Alter default
if r == 0
    r = min(size(im0))/5;
end

% Apply mask
%im0 = applyMask(im0,iC.xTank,iC.yTank);

% Set starting roi
roi = [size(im0,2)/2-r(1)/2 size(im0,1)/2-r(1)/2 2*r 2*r];


% Display options
disp(' ')
disp('COMMANDS  ------------------------------------------------ ')
disp('    left click  : select position in current frame')
disp('    delete      : delete current position')
disp('    right arrow : advance video')
disp('    left arrow  : go back in video')
disp('    1 - 9       : jump to relative position in video segment')
disp('    +           : zoom in')
disp('    -           : zoom out')
disp('    q           : quit interaction mode')
disp(' ')
 
% Create figure window
[hFig, hAxes] = createFigureAndAxes;

% Add data to figure
[hFig, hAxes] = putData(hFig, hAxes);

% Wait for completion of interactive mode
waitfor(hFig)

if isempty(H)
    return
else
    x = H.x;
    y = H.y;
end

    function [hFig, hAxes] = createFigureAndAxes
        
        
        % Close figure opened by last run
        figTag = 'CVST_VideoOnAxis_9804532';
        close(findobj('tag',figTag));
        
        AR = size(im0,2)/size(im0,1);
        %AR = range(iC.yTank)/range(iC.xTank);
        AR_roi = roi(4)/roi(3);
        
        % Relative width of figure window
        rel_width = 0.4;
        
        % Create new figure
        hFig = figure('numbertitle', 'off', ...
            'name', 'Analysis GUI', ...
            'menubar','none', ...
            'toolbar','none', ...
            'resize', 'on', ...
            'tag',figTag, ...
            'renderer','painters',...
            'Units','pixels');
        
        scrsize = get(0,'screensize');
        set(hFig,'Position',scrsize);
        
        %set(hFig,'Units','normalized');
        
        % Size of width and height of 2 windows
        win_size = [scrsize(3)*rel_width scrsize(3)*rel_width/AR];
        
        % Spacing between windows
        spacer    = 10;
        
        % Relative position of main panel
        Lrect(3) = win_size(1);
        Lrect(4) = win_size(2);
        Lrect(1) = scrsize(3)/2 - Lrect(3) - spacer/2;
        Lrect(2) = scrsize(4)/2 - Lrect(4)/2;
        
        Ltitle(3) = Lrect(3);
        Ltitle(4) = 20;
        Ltitle(1) = Lrect(1);
        Ltitle(2) = Lrect(2) + Lrect(4);    
        
        % Relative position of small panel
        Rrect(3) = win_size(1);
        Rrect(4) = win_size(1)/AR_roi;
        Rrect(1) = scrsize(3)/2 + spacer/2;
        Rrect(2) = scrsize(4)/2 - Rrect(4)/2;
        
        Rtitle(3) = Rrect(3);
        Rtitle(4) = 20;
        Rtitle(1) = Rrect(1);
        Rtitle(2) = Rrect(2) + Rrect(4) ;
        
        % Get figure position
        %fPos = get(hFig,'Position');
        
        % Create axes and titles
        hAxes.axis1 = createPanelAxisTitle(hFig,...
            Lrect, Ltitle,'Full frame', 'title1'); % [X Y W H]
        
        hAxes.axis2 = createPanelAxisTitle(hFig, ...
            Rrect, Rtitle, 'ROI', 'title2');

        % Keystrokes
        set(hFig, 'WindowButtonDownFcn', {@butDown, hFig, hAxes});
        
        % Set callback for button press
        set(hFig, 'WindowKeyPressFcn', {@keyPress, hFig, hAxes});

    end

%% Deposite data into figure

    function [hFig, hAxes] = putData(hFig, hAxes)
        
        % Store inital values into axis1
        H.roi = roi;
        H.iFrame = 1;
        H.frames = frames;
        H.vid_path = vid_path;
        H.v = v;
        H.clr = mClr;
        H.x = xStart;
        H.y = yStart;
        H.imInvert = imInvert;

        % Store coordinate data
        guidata(hAxes.axis1, H);
        
        % Render images in figure
        update_fig(hFig, hAxes);
    end


%% Create Axis and Title
% Axis is created on uipanel container object. This allows more control
% over the layout of the GUI. Video title is created using uicontrol.
    function hAxis = createPanelAxisTitle(hFig, posWin, posTitle, axisTitle, textTag)

        % Create panel
        %hPanel = uipanel('parent',hFig,'Position',pos,'Units','Normalized');
        hPanel = uipanel('parent',hFig,'Units','pixels','Visible','on');
        
        % Set position
        hPanel.Position = posWin;
        hPanel.Units = 'pixels';
        
        % Create axis   
        hAxis = axes('Parent',hPanel,'Units','pixels'); 
        hAxis.Position = [0 0 posWin(3:4)];
        %hAxis.Position = [-min(iC.xTank)/2 min(iC.yTank)/2 range(iC.xTank) range(iC.yTank)];
        
        hAxis.XTick = [];
        hAxis.YTick = [];
        hAxis.XColor = [1 1 1];
        hAxis.YColor = [1 1 1];
        
        % Revert units to pixels
        hAxis.Units = 'pixels';
        
        % Set video title using uicontrol. uicontrol is used so that text
        % can be positioned in the context of the figure, not the axis.
        titlePos = posTitle;
        hUI = uicontrol('style','text',...
            'String', axisTitle,...
            'Units','pixels',...
            'Parent',hFig,...
            'Position', titlePos,...
            'BackgroundColor',hFig.Color, ...
            'HorizontalAlignment','left',...
            'Tag',textTag);
        
        % Revert units to pixels
        %hUI.Units = 'pixels';
    end

end



%% Small helper functions

function showFrameOnAxis(hAxis, frame, zoomlevel)
    % This helper function  displays a frame of video on a user-defined axis.

    frame = convertToUint8RGB(frame);

    try
        hChild = get(hAxis, 'Children');
    catch %#ok<CTCH>
        return; % hAxis does not exist; nothing to draw
    end

    isFirstTime = isempty(hChild);

    if isFirstTime
        hIm = displayImage(hAxis, frame);
        zoom(hAxis,zoomlevel)
        if 0
            addScrollPanel(hAxis, hIm);
        end
    else
        hIm = hChild(end);

        try
            set(hIm,'cdata',frame); 
            drawnow;
        catch  %#ok<CTCH>
            % figure closed
            return;
        end
    end
end

function frame = convertToUint8RGB(frame)
    % Convert input data type to uint8
    if ~isa(class(frame), 'uint8')
        frame = im2uint8(frame);
    end

    % If the input is grayscale, turn it into an RGB image
    if (size(frame,3) ~= 3) % must be 2d
        frame = cat(3,frame, frame, frame);
    end
end

function hIm = displayImage(hAxis, frame)
% Display image in the specified axis
frameSize = size(frame);
xdata = [1 frameSize(2)];
ydata = [1 frameSize(1)];
cdata = frame;
cdatamapping = 'direct';

hIm = image(xdata,ydata,cdata, ...
           'BusyAction', 'cancel', ...
           'Parent', hAxis, ...
           'CDataMapping', cdatamapping, ...
           'Interruptible', 'off');
set(hAxis, ...
    'YDir','reverse',...
    'TickDir', 'out', ...
    'XGrid', 'off', ...
    'YGrid', 'off', ...
    'PlotBoxAspectRatioMode', 'auto', ...
    'Visible', 'off');
end

% function [x, y] = full_to_roi(xVal,yVal,roi)
% % Converts from full frame coordinates to ROI coords    
% 
%      % Normalize to size of rendered frame
%      xNorm = (xVal-rect(1))/roi(3);
%      yNorm = (yVal-rect(2))/roi(4);
% 
%      % Transform into frame coords
%      x = xNorm .* range(hAxes.axis2.XLim) + hAxes.axis2.XLim(1);
%      y = yNorm .* range(hAxes.axis2.YLim) + hAxes.axis2.YLim(1);
% end


%% Key press callback

function keyPress(fig, key, hFig, hAxes)
            
    closefig = 0;

     % Load data       
     H = guidata(hAxes.axis2);
     
     % QUIT ('q')
      if strcmp(key.Key,'Q') || strcmp(key.Key,'q')
         
         assignin('caller','H',H);
         
         closefig = 1;
         
     % DELETE
     elseif strcmp(key.Key,'backspace')
         
       % Delete data
       H.x(H.iFrame) = nan;
       H.y(H.iFrame) = nan;
       
        % Store coordinate data
        guidata(hAxes.axis1, H);
       
       % Update figure
       update_fig(hFig, hAxes)
     
     % PLUS (zoom in)
     elseif strcmp(key.Key,'hyphen') || strcmp(key.Key,'minus')
         
         xCntr = H.roi(1) + H.roi(3)/2;
         yCntr = H.roi(2) + H.roi(4)/2; 
         
         % Enlarge roi
         H.roi(3) = 1.5*H.roi(3);
         H.roi(4) = 1.5*H.roi(4);
         H.roi(1) = xCntr - H.roi(3)/2;
         H.roi(2) = yCntr - H.roi(4)/2;
         
        % Store coordinate data
        guidata(hAxes.axis1, H);
       
       % Update figure
       update_fig(hFig, hAxes)
       
     % PLUS (zoom out)
      elseif strcmp(key.Key,'equal') || strcmp(key.Key,'plus')
          
         
         xCntr = H.roi(1) + H.roi(3)/2;
         yCntr = H.roi(2) + H.roi(4)/2; 
         
         % Enlarge roi
         H.roi(3) = H.roi(3)/1.5;
         H.roi(4) = H.roi(4)/1.5;
         H.roi(1) = xCntr - H.roi(3)/2;
         H.roi(2) = yCntr - H.roi(4)/2;
         
        % Store coordinate data
        guidata(hAxes.axis1, H);
       
       % Update figure
       update_fig(hFig, hAxes)
         
     % The following involve changing frame number
     else
         newIndex = [];
         
         % Index of current frame
%          iFrame = find(H.cFrame==H.frames,1,'first');
         
         % RIGHT ARROW/SPACE
         if strcmp(key.Key,'rightarrow') || strcmp(key.Key,'space')
             
             if H.iFrame==length(H.frames)
                 beep
             else
                 % Advance frame
                 newIndex = H.iFrame+1;
             end
             
         % LEFT ARROW
         elseif strcmp(key.Key,'leftarrow')
             
             if H.iFrame==1
                 beep
             else
                 % Reverse frame
                 newIndex = H.iFrame-1;
             end
             
         % NUMBER
         elseif length(key.Key)==1 && (sum(key.Key==num2str([1:9]))==1)
             
             % Requested relative number
             req_num = str2num(key.Key);
             
             % Start 
             if req_num == 1
                 newIndex = 1;
                 
             % End
             elseif req_num==9
                 newIndex = length(H.frames);
                 
             % Something between
             else
                 % Set new frame
                 newIndex = round((req_num/10)*length(H.frames));
             end
                     
         end

         % If new newframe
         if ~isempty(newIndex)

             % Update frame number
             H.iFrame = newIndex;
             
             % Store coordinate data
             guidata(hAxes.axis1, H);
         end
         
      end
     
    if closefig==1
         close(hFig)
    else
      % Update figure
       update_fig(hFig, hAxes)
    end
end


%% Button down callback
function butDown(fig, key, hFig, hAxes)
              
    % Load data       
    H = guidata(hAxes.axis2);   
     
    % Collect current coordinate
    C2 = get(hAxes.axis2, 'CurrentPoint');
    C1 = get(hAxes.axis1, 'CurrentPoint');
    
    hold on
     % If in full frame box . . .
    if (C1(1,1)>=0) && (C1(1,1)<=hAxes.axis1.XLim(2)) && ...
       (C1(1,2)>=0) && (C1(1,2)<=hAxes.axis1.YLim(2))     
        
       % Update roi
       H.roi(1) = C1(1,1)-H.roi(3)/2;
       H.roi(2) = C1(1,2)-H.roi(4)/2;
                
    % If in roi box . . .
    elseif (C2(1,1)>=0) && (C2(1,1)<=hAxes.axis2.XLim(2)) && ...
           (C2(1,2)>=0) && (C2(1,2)<=hAxes.axis2.YLim(2))     
    
       % Current coordinate
       xCurr = H.roi(1) + C2(1,1);
       yCurr = H.roi(2) + C2(1,2);
       
       % Store data
       H.x(H.iFrame,1) = xCurr;
       H.y(H.iFrame,1) = yCurr;
       
    else
        set(gcf,'Pointer','arrow')
    end
    
    % Store coordinate data
    guidata(hAxes.axis1, H);
    
    % Update figure
    update_fig(hFig, hAxes);
             
end


%% Update GUI
    function update_fig(hFig, hAxes)
        
        % Activate figure window
        figure(hFig)
        
        % Load data
        H = guidata(hAxes.axis1);
        
        % List frame number in title 1
        obj = findobj('tag','title1');
        t_str = ['Full Frame :  Frame ' num2str(H.frames(H.iFrame)) ];
        set(obj,'String',t_str);
        
        % List mode in title 2
        obj = findobj('tag','title2');
        set(obj,'String','ROI');
        
        % Index for current frame
%         iFrame = find(H.cFrame==H.frames,1,'first');
      
        % Current point
        xCurr = H.x(H.iFrame);
        yCurr = H.y(H.iFrame);

        % All points
        xAll = H.x(~isnan(H.x));
        yAll = H.y(~isnan(H.y));
        
        % Current frame
        cFrame = H.frames(H.iFrame);
        
        % Read input video frame
        frame = getFrame(H.vid_path,H.v,cFrame,H.imInvert,'gray');

        % Display full video frame
        delete(hAxes.axis1.Children)
        showFrameOnAxis(hAxes.axis1, frame, 0);
        
%         % Get window dimensions (full frame FOR)
%         winWidth1  = size(hAxes.axis1.Children(end).CData,2);
%         winHeight1 = size(hAxes.axis1.Children(end).CData,1);
%         
%         
        % Get coordinates for the ROI (Full frame FOR)
        %rect = get_val;
        
        % Center of ROI (Full frame FOR)
        xCntr = H.roi(1) + H.roi(3)/2;
        yCntr = H.roi(2) + H.roi(4)/2;

        
        % Hold on Axis1
        set(hAxes.axis1,'NextPlot','Add')
        
        % Set H.roi square on full frame
        x1 = H.roi(1); x2 = H.roi(1)+H.roi(3);
        y1 = H.roi(2); y2 = H.roi(2)+H.roi(4);
        plot([x1 x2 x2 x1 x1], [y1 y1 y2 y2 y1], 'w-',...
            mean([x1 x2]).*[1 1],[y1 y2], 'w:', ...
            [x1 x2], mean([y1 y2]).*[1 1], 'w:',...
            'Parent',hAxes.axis1);
        
        % Hold off Axis1
        set(hAxes.axis1,'NextPlot','Replace');
        
        % Get cropped image
        im2 = imcrop(frame, H.roi);
        
        % Display cropped image
        delete(hAxes.axis2.Children);
        showFrameOnAxis(hAxes.axis2, im2, 0);
        
        % Limits for x and y axes
        xL = hAxes.axis2.XLim;
        yL = hAxes.axis2.YLim;
        
        % Hold on Axis 2
        set(hAxes.axis2,'NextPlot','Add');
        
        % Add cross hairs
%         h = plot([mean(xL) mean(xL)],yL,'w:',xL,[mean(yL) mean(yL)],'w:',...
%             'Parent',hAxes.axis2);

        h2 = line(xAll-H.roi(1),yAll-H.roi(2),...
                  'Parent',hAxes.axis2,'Color',[H.clr 0.3],'LineWidth',3);
            
        h2 = scatter(xCurr-H.roi(1),yCurr-H.roi(2),'Parent',hAxes.axis2,...
            'MarkerEdgeColor',H.clr,'MarkerfaceColor','none','SizeData',150);  
        
%         h2 = plot(xCurr-H.roi(1),yCurr-H.roi(2),'+',...
%                   'Parent',hAxes.axis2,'Color',H.clr);
        
        
        % Hold off Axis 2
        set(hAxes.axis2,'NextPlot','Replace');
        
        % Activate figure
        figure(hFig)
        
        % Activate ROI panel
        set(hAxes.axis2,'Selected','on')

    end