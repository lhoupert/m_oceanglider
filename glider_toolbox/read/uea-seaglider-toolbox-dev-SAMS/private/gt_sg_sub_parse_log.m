function data_log = gt_sg_sub_parse_log(glidernum,direct,dives)
% 
% data_log = gt_sg_sub_parse_log(glidernum,direct,dives)
% Reads data from .log files and concatenates them into a single structure
% for integration into the GT_SG toolbox output.
% Is called by gt_sg_load_merge_data.
%
% Inputs:
% glidernum = Glider number
% direct = Working dorectory where .eng and .log files reside
% dives = Dives to process through the toolbox
%          (default: all available)
%
% Outputs:
% data_log = log structure array containing all data in element by element organization
%
% B.Y.QUESTE Feb 2015

% List files in directory
lcfils = dir([direct, 'p', num2str(glidernum), '*.log']);
if isempty(lcfils)
    gt_sg_sub_echo({['ERROR: No .log files found in ',direct,'.'],'Please check directory.'});
    data_log = NaN;
    return
end

% Pull file numbers
dive_nums = cellfun(@(x) str2double(x(5:8)), {lcfils.name});

% Remove launch routine sim dive (dive 0000) and
% reorder structure indexing to account for missing files.
files(dive_nums(dive_nums ~= 0)) = lcfils(dive_nums ~= 0);

% Look for requested dives
if isempty(dives)
    dives = dive_nums(dive_nums ~= 0);
elseif ~isempty(intersect(setxor(dives,dive_nums),dives))
    gt_sg_sub_echo({'ERROR: Requested dives not found.','Missing .log files:',num2str(intersect(setxor(dives,dive_nums),dives))});
    data_log = NaN;
    return
end

% Initialise waitbar
gt_sg_sub_echo(['Parsing .log files in ',direct]);
h = waitbar(0,'Parsing .log files...') ;
count = numel(dives);

% Initialise log structure array
data_log(max(dive_nums)).glider = [];

% Pull data
for dstep = dives
    data_log = gt_sg_sub_append_to_struct(data_log,getdata([direct files(dstep).name]),dstep);
    waitbar(dstep/count,h)
end
close(h) % Close waitbar
fclose('all');
end

function data = getdata(filename)

fid = fopen(filename,'r');

data = [];
if fid ~= -1
    
    % Parse version, glider number and start time etc.
    linedata = fgetl(fid);
    line_num = 1;
    
    while ~any(strfind(linedata,'data:')) && ischar(linedata)
        value = linedata(strfind(linedata,':')+1:end);
        % Don't use str2double as it converts strings with spaces to NaNs
        % instead of arrays.
        data.(linedata(1:strfind(linedata,':')-1)) = str2num(value(regexp(value,'[0-9, ]')));
        linedata = fgetl(fid);
        line_num = line_num+1;
    end
    
    % Skip line beginning "data:"
    linedata = fgetl(fid);
    line_num = line_num+1;
    
    % Parse each field sequentially
    while linedata ~= -1
        try
            delimiter = strfind(linedata,',');
            fieldname = linedata(2:delimiter-1);
            
            % Computationally more expensive to test if fieldname needs to 
            % be changed with genvarname than parsing each one by default.
            % So leaving it as it is for now...
            fieldname(strfind(fieldname,'%')) = [];
            fieldname(strfind(fieldname(1),'_')) = [];
            fieldname(strfind(fieldname(1),'$')) = [];
            if regexp(fieldname(1),'[0-9]')
                fieldname = ['x' fieldname];
            end
            
            fdata = str2num(linedata(delimiter+1:end));
            
            % Initialise field for future concatenation
            if ~isfield(data,fieldname)
                data.(fieldname) = [];
            end
            
            if ~isempty(fdata) % If data is num (ie. not failed str2num), store it.
                data.(fieldname) = [data.(fieldname); fdata;];
            else % If data is a string (ie. failed str2num), put string instead.
                data.(fieldname) = [data.(fieldname); {linedata(strfind(linedata,',')+1:end)};];
            end
        catch
            gt_sg_sub_echo({['ERROR: Could not parse ',filename,' at line ',num2str(line_num),'.'],['Please verify the integrity of the .log file'],['"',linedata,'"']});
        end
        linedata = fgetl(fid);
        line_num = line_num+1;
    end
else
    data = NaN;
end

fclose(fid);

end