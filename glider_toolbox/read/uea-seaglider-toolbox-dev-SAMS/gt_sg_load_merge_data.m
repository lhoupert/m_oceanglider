function data = gt_sg_load_merge_data(glidernum,direct,dives,sg_cal_direct)
%
% data = gt_sg_load_merge_data(glidernum,direct,start_profile,count,sg_cal_direct)
%
% Load Seaglider log and engineering data into the GT_SG toolbox data format 
% for further processing.
%
% Inputs:
% glidernum = Glider number
% direct = Working dorectory where .eng and .log files reside
%          (default: current directory)
% dives = Dives to process through the toolbox
%          (default: all available)
% sg_cal_direct = Directory of sg_calib_constants.m
%          (default: current directory)
%
% Outputs:
% data = GT_SG structure containing all data in element by element organization
%
% B.Y.QUESTE Feb 2015

% Parse inputs
% Require glider number
if isempty(nargin)
    disp('Please provide a Glider number.');
    data = NaN;
    return
end

% Specifiy working directory
if nargin < 2
    direct = [pwd filesep];
end
if ~strcmp(direct(end),filesep)
    direct = [direct filesep];
end

% Specify sg_calib_constants directory
if nargin < 4
    sg_cal_direct = direct;
end
if ~strcmp(sg_cal_direct(end),filesep)
    sg_cal_direct = [sg_cal_direct filesep];
end

% Initialise output data variable
data = [];
gt_sg_sub_echo(['BEGINNING PROCESSING OF SG',num2str(glidernum)]);

% Load sg_calib_constants
if exist([sg_cal_direct,'sg_calib_constants.m'],'file')
    data.sg_calib_const = gt_sg_import_constants([sg_cal_direct 'sg_calib_constants.m']);
else
    gt_sg_sub_echo({'ERROR: No sg_calib_constants.m found',...
    'Please provide path to the file',...
    'or add it to the current directory.',...
    '.log and .eng files will be still be parsed.'});
end

% Load data
if ~exist('dives','var')
    dives = [];
end

data0.log = gt_sg_sub_parse_log(glidernum,direct,dives);
data0.eng = gt_sg_sub_parse_eng(glidernum,direct,dives);

%----------------------------------------------------------------------------------------------------------------------
% Remove data with nan for condFreq or tempFreq or negative depth before detecting for single level dive % L Houpert 25/01/201
data0 = gt_sg_sub_remove_sbect_nan(data0,1);
% remove profile with only 1 data point (solve problems in gt_sg_sub_diff. Work for sg605-m1)
ccbad= [];
for ijk=1:length(data0.eng)
    if length(data0.eng(ijk).elaps_t)<=1 & ~isempty(data0.eng(ijk).columns) 
        ccbad = [ccbad data0.eng(ijk).dive];
        disp(['Not enough data points for dive : '  num2str(data0.eng(ijk).dive)])
        disp(['Dive removed from the list of the dives to process'])     
    end
end
if length(ccbad)>0
	data.log = gt_sg_sub_parse_log(glidernum,direct,setdiff(dives,ccbad));
	data.eng = gt_sg_sub_parse_eng(glidernum,direct,setdiff(dives,ccbad));
else
	data.log = data0.log;
	data.eng = data0.eng;	
end

%--------------------------------------------------------------------------------------
% After 2nd parsing, remove data with nan for condFreq or tempFreq or negative depth % L Houpert 25/01/2016
data = gt_sg_sub_remove_sbect_nan(data);

if ~isstruct(data.eng) && ~isstruct(data.log)
    gt_sg_sub_echo({['ERROR: no dives found in ',direct,'.'],'Aborting.'});
    return;
end

% Diplay output metrics
missing_pairs = setxor([data.eng.dive],[data.log.dive]);

if numel(missing_pairs) > 0
    gt_sg_sub_echo({'WARNING: Missing matching .eng and .log files.',['Missing .log files for dives: ',num2str(intersect(missing_pairs,[data.eng.dive]))],['Missing .eng files for dives: ',num2str(intersect(missing_pairs,[data.log.dive]))]});
end

gt_sg_sub_echo({['Processed ',num2str(sum(~isnan([data.log.dive]))),' .log and ',num2str(sum(~isnan([data.eng.dive]))),' .eng files succesfully.']});

% Output as columns (because it's prettier)
data.log = data.log(:);
data.eng = data.eng(:);
end
