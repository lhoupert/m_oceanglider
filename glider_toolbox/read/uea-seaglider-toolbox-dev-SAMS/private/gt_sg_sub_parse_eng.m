function data_eng = gt_sg_sub_parse_eng(glidernum,direct,dives)
% 
% data_eng = gt_sg_sub_parse_eng(glidernum,direct,dives)
% Reads data from .eng files and concatenates them into a single structure
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
% data_eng = eng structure array containing all data in element by element organization
%
% B.Y.QUESTE Feb 2015

% List files in directory
lcfils = dir([direct, 'p', num2str(glidernum), '*.eng']);
if isempty(lcfils)
    gt_sg_sub_echo({['ERROR: No .eng files found in ',direct,'.'],'Please check directory.'});
    data_eng = NaN;
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
    gt_sg_sub_echo({'ERROR: Requested dives not found.','Missing .eng files:',num2str(intersect(setxor(dives,dive_nums),dives))});
    data_eng = NaN;
    return
end

% Initialise waitbar
gt_sg_sub_echo(['Parsing .eng files in ',direct]);
h = waitbar(0,'Parsing .eng files...') ;
count = numel(dives);

% Initialise eng structure array
data_eng(max(dive_nums)).glider = [];

% Pull data
for dstep = dives
    data_eng = gt_sg_sub_append_to_struct(data_eng,getdata([direct files(dstep).name]),dstep);
    waitbar(dstep/count,h)
end
close(h) % Close waitbar

% Apply same substructure format to empty dives
% Format copied from Sunke, so don't judge me for the sloppiness.
emptyDives = setxor([1:numel(data_eng)],dives);
if ~isempty(emptyDives)
    fnames = fieldnames(data_eng);
    for fstep = 1:numel(fnames)
        
        cellArray = {data_eng.(fnames{fstep})};
        isSubstructure = cellfun(@(x) isstruct(x),cellArray);
        if any(isSubstructure)
            substructureInd = find(isSubstructure,1,'first');
            fnamesSub = fieldnames(data_eng(substructureInd).(fnames{fstep}));
            for fstep_sub = 1:numel(fnamesSub)
                for istep = find(~isSubstructure)
                    data_eng(istep).(fnames{fstep}).(fnamesSub{fstep_sub}) = [];
                end
            end
        end
    end
end

fclose('all');
end


function data = getdata(filename)

% Open file
fid = fopen(filename,'r');

% Return empty arrray if fail
data = [];

if fid ~= -1 % If file exists
    
    % Parse version, glider number and start time etc.
    linedata = fgetl(fid);
    line_num = 1;
    while ~any(strfind(linedata,'data:')) && ischar(linedata)
        data.(linedata(2:strfind(linedata,':')-1)) = strtrim(linedata(strfind(linedata,':')+1:end));
        linedata = fgetl(fid);
        line_num = line_num+1;
    end
    
    % Convert dive number from string to double
    data.dive = str2double(data.dive);
    
    % Convert start date from string to array
    data.start = str2num(data.start);
    
    % Parse data column header
    data.columns = strsplit(data.columns,',');
    
    % Import data
    dataArray = textscan(fid, repmat('%f',1,numel(data.columns)), 'Delimiter', ' ', 'ReturnOnError', false);
    
    % Ok, ok, I know eval is sloppy, but it's pretty bloody quick this way,
    % and it copes really well with the multi-level substructures...
    % Transpose so multi-dive concatenatable with [].
    for istep = 1:numel(dataArray)
        eval(['data.',data.columns{istep},' = dataArray{istep}'';']);
    end 
end

fclose(fid);
end
