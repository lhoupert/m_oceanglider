function [dataf]=sig_filt(window,data,wsize)
%[dataf]=sig_filt(window,data,wsize);
if rem(wsize,2)~=0; wsize=wsize+1;end
    
Mdata=filter(window,1,data);

ifen=floor(wsize/2);

dataf=[data(1:(wsize+ifen),:); Mdata(2*wsize+1:end-(2*wsize+1),:);...
    data(end-(2*wsize+ifen):end,:)];

end
