function dx = gt_sg_sub_diff(x)
%
% dx = gt_sg_sub_diff(x)
% First order differentiation returning an array of the same length.
% In an array is input, the it will return dx, if a cell array containing
% two same length arrays is input, it will return dx{1}/dx{2}.
%
% Inputs:
% x = array or 2x1 cell array containing two arrays of the same length.
%
% Outputs:
% dx = first order differential of x
%
% B.Y.QUESTE Feb 2015
% S. SCHMIDTKO 2013
%
if iscell(x)
    
    in1 = x{1};
    in2 = x{2};
    
    % TODO: This is Sunke's method. Verify it makes sense. Produces
    % different results to diff etc (slightly smoothed). Is this an issue?
    d1 = convn(in1,[1 0 -1]./2,'same'); % central differences
    d1([1 end]) = [diff(in1([1 2])) diff(in1([end-1 end]))]; % fix end members
    
    d2 = convn(in2,[1 0 -1]./2,'same'); % central differences
    d2([1 end]) = [diff(in2([1 2])) diff(in2([end-1 end]))]; % fix end members
    
    dx = d1./d2;

else
    
    dx = convn(x,[1 0 -1]./2,'same'); % central differences
    dx([1 end]) = [diff(x([1 2])) diff(x([end-1 end]))]; % fix end members
    
end
end