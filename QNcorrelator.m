function [Act, R, D, S] = QNcorrelator(phi_cl,qgrid,xi,active,ell);
% Negative action of Q(x_i), i=1..N.  refer to Bialek et.al. 1996,
% Eq.(19) 
% IN: 
%   phi_cl -- calculated classical solution 
%   qgrid -- grid 
%   xi -- samples 
%   ell -- smoothness scale 
% OUT: 
%   Act -- negative action
%   R -- fluctuation determinant
%   D -- prior term
%   S -- data term

% the incremental step
dx=qgrid(2)-qgrid(1);
% number of sample points
N=length(xi);
scale = N/sum(active);

% spline for phiclass
disp('   Creating spline representation for the current classical solution');
phispline = spline([qgrid,qgrid(end)+dx], [phi_cl,phi_cl(1)]);

% part of action coming from data samples
disp('   Evaluating classical source term.');
S= sum(ppval(phispline, xi));
% part of action coming due to smoothness
disp('   Evaluating classical prior term.');
D= ell/2*dx*sum((cdiffl(phi_cl)/dx).^2);
% fucntional determinant
disp('   Evaluating functional determinant');
R=Rfunctdet(phi_cl, active, scale, dx,ell);

% the correlator (action)
Act=S+D-R;