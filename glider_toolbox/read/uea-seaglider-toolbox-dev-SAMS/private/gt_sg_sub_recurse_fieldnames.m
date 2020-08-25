function nameList = gt_sg_sub_recurse_fieldnames(data,varargin)

if nargin < 3
    nameList = {};
else
    nameList = varargin{2};
end
    
if nargin < 2
    dataName = inputname(1);
else
    dataName = varargin{1};
end
    
if isstruct(data)
    fieldNames = fieldnames(data(1));
    nFields = numel(fieldNames);
    for jstep = 1:nFields
        nameList = gt_sg_sub_recurse_fieldnames(data(1).(fieldNames{jstep}),[dataName,'.',fieldNames{jstep}],nameList);
    end
else
    nameList{end+1} = dataName;
    return;
end
end