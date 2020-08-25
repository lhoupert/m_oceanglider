% ********************************************************************************
%
% ********************************************************************************
function x = fillchk(struct, field, len)
    fn = fieldnames(struct);

    for ii = 1:length(fn)
        fname = fn{ii};
        sr_new.(lower(fname)) = struct.(fname);
    end 
 
    if isfield(sr_new, lower(field))
        x = getfield(sr_new, char(lower(field)));
    else
        % RETURN AN EMPTY SET SO THAT WE KNOW THAT NO DATA EXISTS!
        x = [];
    end
