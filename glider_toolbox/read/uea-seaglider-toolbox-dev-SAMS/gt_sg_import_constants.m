function constants = gt_sg_import_constants(file)
%
% calibs = gt_sg_import_constants(direct)
%
% Loads data from .m files (ie. sg_calib_constants) to GT_SG toolbox format.
%
% Inputs:
% direct = Path to sg_calib_constants.m (or other)
%
% B.Y.QUESTE Feb 2015

% Throw error if file cannot be found.
if ~exist(file,'file')
    gt_sg_sub_echo({['ERROR: ',file,' not found.']});
    constants = [];
    return;
end

% Parse file name and path
[path_to_constants,file_name,file_ext] = fileparts(file);

% Find file in matlab seach path if path not provided (ie. default)
if isempty(path_to_constants)
    path_to_constants = fileparts(which(file));
end

% Go to file for easy loading
cpath = cd(path_to_constants);

% Decline all responsibilities - let the user make sure it's the right file
gt_sg_sub_echo(['Reading ',file_name,file_ext,' from directory: ' path_to_constants]);

% Load sg_calib_constants to this workspace
eval(file_name)

% Create a full path to save in constants for future reference/verification
file_path = [path_to_constants,filesep,file_name,file_ext];

% Return to original working directory
cd(cpath)

% Clear variables we don't want to save in constants sub-structure
clear cpath path_to_constants file_name file_ext file

% Save to then load from temporary matfile to move to structure format.
tmpname =  [tempname '.mat'];
save(tmpname);
constants = load(tmpname);
constants = rmfield(constants,'tmpname');
delete(tmpname);

end

