function out = gt_sg_sub_filter(in,range_median,range_lowpass)
%
% x = gt_sg_sub_filter(x)
% Running median filter with window of "2*range_median + 1" range
% followed by a lowpass box filter of "2*range_median + 1" range.
%
% Inputs:
% in = array of data
%
% Outputs:
% out = Filtered and smoothed x
%
% B.Y.QUESTE Feb 2015
% S. SCHMIDTKO 2013
%

s = size(in);
in = in(:);
n = numel(in);
m = nan(n+range_median*2,range_median*2+1);

m(range_median+1:end-range_median,:) = repmat(in,1,range_median*2+1);

for istep=1:range_median*2+1
    m(istep:istep+n-1,istep) = in;
end
out = nanmedian(m,2);

out = out(range_median+1:end-range_median);

out = conv(out,ones(1,range_lowpass*2 +1)/(range_lowpass*2 +1),'same');

if size(out) ~= s
    out = out';
end

end