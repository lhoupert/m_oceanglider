function [dataf]=sig_filt_lp(Fc,data,Nn);
% [dataf]=sig_filt_lp(Fc,data,Nn);
%  input: data: data to be filtered
%         Nn:   butterworth filter order 
%         Fc:   cutoff frequency

if rem(Nn,2)~=0; Nn=Nn+1;end
ifen=floor(Nn/2);

Hf = fdesign.lowpass('N,Fc',Nn,Fc);
Hd(1) = design(Hf,'window','window',@butter);
Mdata=filter(Hd(1),data);
% c=0;
% while c<4
%     Mdata=[data(1:(Nn+ifen),:); Mdata(2*Nn+1:end-(2*Nn+1),:);...
%     data(end-(2*Nn+ifen):end,:)];
%     Mdata=filter(Hd(1),Mdata);   
%     c=c+1;    
% end



dataf=[data(1:(Nn+ifen),:); Mdata(2*Nn+1:end-(2*Nn+1),:);...
    data(end-(2*Nn+ifen):end,:)];

end
