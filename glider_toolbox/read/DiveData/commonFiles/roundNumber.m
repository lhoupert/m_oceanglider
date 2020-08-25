% *************************************************************************
% Purpose: the purpose of this function is to round numbers to the
%          specified decimal digits of precision
% *************************************************************************
function [nRound]= roundNumber(nNumber, nPrecision)
  mult = 10^nPrecision;
  nRound = round(nNumber*mult)./mult;
end
