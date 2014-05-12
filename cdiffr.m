function [out]=cdiffr(x);
%cyclic difference - to right
out=zeros(size(x));
out(1:end-1)=diff(x);
out(end)=x(1)-x(end);

