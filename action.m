function [out]=action(log_ell, xi, xi_hist, active, qgrid, d2, L, opt);
% This function calculates the action for minimization over l's
%   log_ell -- the investigated value of log(l_smooth)
%   xi -- data
%   xi_hist -- histogrammed data
%   active -- activity sequence
%   qgrid -- grid
%   d2 -- d^2 differetial operator
% opt governs the behavior of the function
% opt=0 -- calculate action
%    =1 -- reset the pre-calculated values for action


% we store already calculated values
global l_calc phi_calc Act_calc S_calc D_calc R_calc

ell=exp(log_ell);

Nqgrid=length(qgrid);
if (opt==1)
  % reset the pre-calculated values
  l_calc = 100;
  phi_calc = zeros(Nqgrid,1);
  Act_calc =0;
  D_calc   =0;
  R_calc   =0;
  S_calc   =0;
  out=0;
elseif opt==0
  % calculate the action
  phi  = zeros(1,Nqgrid);
  Act  = 0;
  R    = Act; S = Act; D = Act;
  N    = length(xi);   
  
  % we find the "best" l (indexed by i) among pre-calculated to supply
  % its phi as the initial solution
  [value i]=min((l_calc - ell).^2);
  
  if value==0
    % already worked at this l
    out = Act_calc(i);
  else
    % did not work at this l, go ahead with calculations
    % calculating classical solutions for N(i) points.
    disp(['  Smoothness scale is l_s=' num2str(ell) ...
          '. Investigating ' num2str(N) ' points.']);
    % solving for classical solution at given l
    phi=qclass(xi_hist, active, qgrid, d2, L, ell, phi_calc(:,i));
    [Act,R,D,S]= QNcorrelator(phi', qgrid, xi, active,ell);
    
    l_calc      = [l_calc, ell];
    phi_calc    = [phi_calc, phi];
    Act_calc    = [Act_calc, Act];
    S_calc      = [S_calc, S];
    D_calc      = [D_calc, D];
    R_calc      = [R_calc, R];

    out= Act;
   end
end;
