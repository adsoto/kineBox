function imOut = applyMask(im,xMask,yMask,clr,bwInvert)
% Applies mask to color or grayscale image

if nargin<4
    clr = [255 255 255];
end

if nargin<5
    bwInvert = 0;
end

if ~isempty(xMask)
    % Color
    if size(im,3)==3
        
        % Get individual channels
        imR = im(:,:,1);
        imG = im(:,:,2);
        imB = im(:,:,3);
        
        if bwInvert==1
            % Make binary mask of tank region
            bwMask = roipoly(imR,xMask,yMask);
        else
            % Make binary mask of tank region
            bwMask = ~roipoly(imR,xMask,yMask);
        end
        
        % Set white mask in each
        imR(bwMask) = clr(1);
        imG(bwMask) = clr(2);
        imB(bwMask) = clr(3);
        
        % Package into image
        imOut(:,:,1) = uint8(imR);
        imOut(:,:,2) = uint8(imG);
        imOut(:,:,3) = uint8(imB);
        
    % Grayscale
    else
        % Set default output
        imOut = im;
        
        if bwInvert==1
            % Make binary mask of tank region
            bwMask = roipoly(im,xMask,yMask);
        else
            % Make binary mask of tank region
            bwMask = ~roipoly(im,xMask,yMask);
        end
        
        % Apply mask
        imOut(bwMask) = 255;
        
    end
    
else
    imOut = im;
end
