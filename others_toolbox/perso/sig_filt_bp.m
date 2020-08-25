function [dataf]=sig_filt_bp(Fc1,Fc2,data,Nn)
% [dataf]=sig_filt_lp(Fc,data,Nn);
% Lowpassdata with an Nnth order hamming window filter with a cutof frequencie Fc
if rem(Nn,2)~=0; Nn=Nn+1;end
    
Hf = fdesign.bandpass('N,Fc1,Fc2',Nn,Fc1,Fc2);
Hd(1) = design(Hf,'window','window',@hamming);
Mdata=filter(Hd(1),data);

ifen=floor(Nn/2);

dataf=[data(1:(Nn+ifen),:); Mdata(2*Nn+1:end-(2*Nn+1),:);...
    data(end-(2*Nn+ifen):end,:)];

end
