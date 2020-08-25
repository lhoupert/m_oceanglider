function y = nancumsum(x)
x(isnan(x))=0;
y=cumsum(x);
