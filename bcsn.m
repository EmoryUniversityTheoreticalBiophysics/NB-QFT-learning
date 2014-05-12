function [ell phi active] = bcsn(xi,L,bnds,Nqgrid);
% function [ell phi active] = bcsn(xi,L,bnds);
%  
% Finds the underlying rate for a Poisson process or the underlying
% PDF using the BCS/NB (hereafter, BCSN) method.
%
% In:
%  xi - poisson event times or sample coordinates
%  L - range of the variables, or the total duration of the
%    sequence
%  bnds -- vector many times 2 -- the boundaries of regions
%    (start/end) when the Poisson process is active, or a sample is
%    allowed; can be ommited, that is, set to [] if the process is always active
%  Nqgrid - may be omitted,  default 30000, number of grid points
%
% Out:
%  ell -- inferred smoothness scale
%  phi -- inverred rate or the log of the probability distribution
%  active -- the fraction of each time bin when the process is  active
%
% I am unsure if the code works correctly for L~=1
%
% (c) Ilya Nemenman, ilya@menem.com, 2000-2005

%preparing data: removing NaN's and data out of [0:L) range
xi=xi(~isnan(xi));
xi=xi(xi>=0);
xi=xi(xi<L);

if(~exist('Nqgrid'))
  Nqgrid=30000;
end
[xi_hist,active,qgrid,d2] =histQ(xi,L,Nqgrid,bnds);

%opts =optimset('Display', 'iter', 'TolX', 1e-2);
opts=[];

global l_calc phi_calc
[log_ell]=fminbnd('action',log(1e-4*L),log(100*L),opts,xi,xi_hist, ...
                           active, qgrid, d2, L,0);   
ell=exp(log_ell);
phi=phi_calc(:,find(l_calc==ell));

