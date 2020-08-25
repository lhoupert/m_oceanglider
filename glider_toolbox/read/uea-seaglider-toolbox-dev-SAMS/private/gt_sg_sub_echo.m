function output = gt_sg_sub_echo(message)
%
% gt_sg_sub_echo(message)
%
% Private subfunction of GT_SG toolbox responsible for printing output to
% screen and to a log file.
%
% Inputs:
% message = cell array of strings to print
%
% B.Y.QUESTE Feb 2015

persistent gt_sg_log

fid = fopen('gt_sg_processing.log','a');

if ~iscell(message)
    message = {message};
end

timestamp = datestr(now);
format = ['\ngt_sg >> %s',repmat('\n  >  %s',1,numel(message)),'\n'];

% To variable
gt_sg_log = [gt_sg_log,sprintf(format,timestamp,message{:})];

% To screen
fprintf(format,timestamp,message{:});

% To log file
fprintf(fid,format,timestamp,message{:});
fclose(fid);

output = gt_sg_log;

end