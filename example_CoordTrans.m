function example_CoordTrans
% Code to test the coordinate transformation code


% Body length
bLen = 5;

% Orientation wrt x-axis of global system
theta = pi/8;

% Coordinates for head
head = [3 3];

% Coordindates for tail 
tail(1,1) = head(1) + bLen*cos(theta);
tail(1,2) = head(2) + bLen*sin(theta);

% Coordinates for fins
rightfin = [4 5];
leftfin = [4.2 3.2];

% Define coordinate system with 
tform = defineSystem2d(head,tail,'x');

% Transform into local coordinates
headL     = transCoord2d(head,tform,'global to local');
tailL     = transCoord2d(tail,tform,'global to local');
rightfinL = transCoord2d(rightfin,tform,'global to local');
leftfinL  = transCoord2d(leftfin,tform,'global to local');

% Transform back into global coordinates
head2     = transCoord2d(headL,tform,'local to global');
tail2     = transCoord2d(tailL,tform,'local to global');
rightfin2 = transCoord2d(rightfinL,tform,'local to global');
leftfin2  = transCoord2d(leftfinL,tform,'local to global');

figure
subplot(1,2,1)
plot(head(1),head(2),'ro',[head(1) tail(1)],[head(2) tail(2)],'r-',...
    [head(1) rightfin(1)],[head(2) rightfin(2)],'r--',...
    [head(1) leftfin(1)],[head(2) leftfin(2)],'r--',...
    headL(1),headL(2),'bo',[headL(1) tailL(1)],[headL(2) tailL(2)],'b-',...
    [headL(1) rightfinL(1)],[headL(2) rightfinL(2)],'b--',...
    [headL(1) leftfinL(1)],[headL(2) leftfinL(2)],'b--')
axis equal
title('Global to local')

subplot(1,2,2)
plot(head2(1),head2(2),'ro',[head2(1) tail2(1)],[head2(2) tail2(2)],'r-',...
    [head2(1) rightfin2(1)],[head2(2) rightfin2(2)],'r--',...
    [head2(1) leftfin2(1)],[head2(2) leftfin2(2)],'r--',...
    headL(1),headL(2),'bo',[headL(1) tailL(1)],[headL(2) tailL(2)],'b-',...
    [headL(1) rightfinL(1)],[headL(2) rightfinL(2)],'b--',...
    [headL(1) leftfinL(1)],[headL(2) leftfinL(2)],'b--')
axis equal
title('Local to global')

