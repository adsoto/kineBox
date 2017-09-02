function varargout = imInteract(im,action)
% Image interactive mode 
%   im       - image 
%   action   - strong requesting type of iteraction 
%
% [tVal,x,y] = imInteract(im,'threshold and selection')
%    returns the chosen threshold and blob position
%
% Developed by McHenryLab, UC Irvine


% Prompt for threshold and blob selection
if strcmp(action,'threshold and selection')
    
    % Plot image
    figure
    imshow(im,'InitialMagnification','fit');
    title('Choose threshold w/arrows, select blob, press return')
    hold on
    
    % Guess threshold
    tVal = graythresh(im);
    
    % Initialize position vectors
    xPos = [];
    yPos = [];
    
    % Loop interaction
    while true
        
        % Make binary
        bw = im2bw(im,tVal);
        
        % Fill holes
        bw = imfill(bw,'holes');
        
        if ~isempty(xPos)
            bw = bwselect(bw,xPos,yPos);
        end
        
        % Trace perimeter
        [y, x] = find(bwperim(bw,8));
        
        % Overlay
        h = plot(x,y,'.g');
        
        % Get input
        [x,y,b] = ginput(1);
        
        % Up arrow
        if b==30
            tVal = min([tVal+0.05 1]);
            
        % Down arrow
        elseif b==31
            tVal = max([tVal-0.05 0]);
            
        % Left click
        elseif b==1
            xPos = x;
            yPos = y;
            
        % esc
        elseif b==27   
            
            % Restart coordinates
            xPos = [];
            yPos = [];
            
        % Return
        elseif isempty(b)
            
            % Check for coordinates
            if isempty(xPos)
                warning('You need to select a blob before pressing return')
            else
                break
            end
        end
        
        % Remove perimeter
        delete(h)       
    end
    
    % Function outputs
    varargout{1} = tVal;
    varargout{2} = xPos;
    varargout{3} = yPos;
    
else
    error('requested action not recognized');
end

