function [dataf]=moy_gl0(data,fen);
% function [dataf]=moy_gl(data,fen);

Mdata=filter(ones(1,fen)/fen,1,data);

ifen=round(fen*0.5);

dataf=[data(1:(fen+ifen))' Mdata(2*fen+1:end-(2*fen+1))'...
    data(end-(2*fen+ifen):end)']';

end
