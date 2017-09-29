function paths = givePaths(rootPath)
% Returns a structure of paths for running code

% If running on Matt's laptop . . .
if nargin<1 && ~ispc && ~isempty(dir('/Users/mmchenry/Documents/Matlab code'))
    % Define root in Matt's folder
    rootPath = '/Users/mmchenry';
    
    paths.vid_root = [rootPath '/Documents/Projects/kineBox/video'];
    
    paths.data_root = [rootPath '/Documents/Projects/kineBox/data_seastars'];
end

%TODO: Add your own default, like done above and modify below to match
% your directory strcuture

% Path to sample videos in Google Drive 
%paths.vid_root = [rootPath filesep 'Google Drive' filesep 'Projects' filesep 'kineBox'];


