function [R]=Rfunctdet(phiclass, active, scale, dx, l_s);
% calculates the log of R, functional determinant, for given
% phiclass --	classical solution
% active   --   activity mode
% scale    --	number of samples (scaled), that is N or
%               N/active(T)
% dx       --   step
% l_s      -- 	smoothness scale
% refer to Bialek et.al. 1996, Eq. (21).

% if all "active" are the same --> no inactive regions -->
% can use WKB approximation
if(all(active==active(1)))
  R=(-1/2*sqrt(scale/l_s)*sum(active.*exp(-0.5*phiclass)));
else                                    % otherwise do fill
                                        % numerical integration
  T=sum(active);
  y=zeros(size(phiclass));
  x=zeros(size(phiclass));
  yp=zeros(size(phiclass));
  xp=zeros(size(phiclass));
  
  yp(1)=1;                              % initial condition
  x(1) =Inf;
  do_y=1;
  for i=2:length(phiclass)
    if((y(i-1)<10)&&do_y)
      y(i) = y(i-1) +yp(i-1)*dx;
      yp(i)= yp(i-1)+y(i-1)*active(i)*exp(-phiclass(i))*scale/l_s;
      x(i) = -log(y(i))*sqrt(l_s/scale);
      xp(i)= -yp(i)/y(i)*sqrt(l_s/scale);
      %    disp(['y=' num2str(y(i)) ' yp=' num2str(yp(i))])
    else
      do_y=0; 
      x(i) = x(i-1) +xp(i-1)*dx;
      xp(i)= xp(i-1)+sqrt(scale/l_s)* ...
             (xp(i-1)^2-active(i)*exp(-phiclass(i))/dx)*dx;
      if(abs(xp(i)-xp(i-1))>10)
        error(['Instability occured in numerical integration of van ' ...
               'Vleck ODE.'])
      end
    end
    %  disp(['x=' num2str(x(i)) ' xp=' num2str(xp(i))])
  end
  R = 1/2*sqrt(scale/l_s)*x(end);
end
