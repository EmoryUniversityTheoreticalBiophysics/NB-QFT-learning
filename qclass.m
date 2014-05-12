function [phi]=qclass(xi_hist, active, qgrid, d2, L, l_smooth, phi_in);
% 
% Solving for classical solution  of the learning problem (as in BCS'96).
% IN:
%  xi_hist -- histogrammed samples
%  active -- bin fraction to enter normalization
%  qgrid -- grid
%  d2 -- d^2 differential operator
%  L -- the box size
%  l_smooth -- smoothness scale
%  phi_in -- initial solution
% 
% OUT:
%  phi -- classical solution


% qgrid size and step
Nqgrid = length(qgrid);
dx=L/Nqgrid;
% # of samples and normalization enforcement
N=sum(xi_hist);
scale = N/sum(active);

% needed for periodic wrapping for banded matrices
global bcsn_overhang

phi=phi_in(:);			        % initial value of phi
E2=1;			                % approximation error
j=0;				        % cycle counter
while max(abs(E2))>1e-5                 % natural scale for phi is
                                        % mean phi=0, deviations of
                                        % order 1
  j=j+1;
  if(floor(j/10)==j/10)  	        % progress indication every 10 cycles
    disp(['     Doing iteration number ' num2str(j) ... 
          '; maximal error is ' num2str(max(abs(E2))) ... 
          '. l_s=' num2str(l_smooth) ', N=' num2str(N) '.']);   
  end;   

  emp = scale*exp(-phi).*active(:);     % scale*exp(-phi)*Xi
  % calculating the current value of error  
  E2= l_smooth*conv([1 -2 1]/dx,[phi(end);phi;phi(1)]);
  E2= reshape(E2(3:end-2),[length(phi),1])+emp -xi_hist;

  % check for numerical instability
  if max(abs(E2))>10000, disp('Method may diverge!'), end;

  % calculating the expected change in phi by Newton-Raphson
  if(size(d2,2)>2)                      % matlab
    dphi=(d2*l_smooth -spdiags(emp,0,Nqgrid,Nqgrid))\E2; 
  else                                  % octave
    dphi=SBSolve(d2*l_smooth - [emp(end-bcsn_overhang+1:end); ...
                        emp;emp(1:bcsn_overhang)]*[1 0], ...
                 [E2(end-bcsn_overhang+1:end);E2;...
                  E2(1:bcsn_overhang)]);
    dphi=dphi(bcsn_overhang+1:(end-bcsn_overhang));
  end
  phi=phi-dphi;
end

