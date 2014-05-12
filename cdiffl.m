function [out]=cdiffl(x);
%cyclic difference - to left
out=zeros(size(x));
out(2:end)=x(2:end)-x(1:end-1);
out(1)=x(1)-x(end);

