function [xi_hist active qgrid d2]=histQ(xi, L, Nqgrid, bnds);
% This function histograms data as requested. 
% IN:
%  xi -- samples
%  L -- the box size
%  Nqgrid -- number of grid points
%  bnds -- long times 2 matrix indicating beginning and end of each
%    interval to be counted in normalization
% 
% OUT:
%  xi_hist -- histogrammed data
%  active -- whether particular bin is to be counted in
%     normalization
%  qgrid -- grid values
%  d2 -- second derivative operator

dx=L/Nqgrid;

% are we in MatLab or Octave?
a=version;
if(length(findstr(a,'R'))>0)
  matlab=1;
else
  matlab=0;
end

% The qgrid is cyclic. That is, qgrid(Nqgrid+1) is
% the same as qgrid(1); that is, L=0, and last qgrid node is <L
qgrid = [0:dx:L];
xi_hist=hist(xi, qgrid);                % histogram of samples, 
xi_hist=[xi_hist(1)+xi_hist(Nqgrid+1),xi_hist(2:Nqgrid)];     
xi_hist=xi_hist(:);                     % this is now the 
                                        % # of samples thet fell in 
                                        % qgrid point+-dx/2
qgrid = qgrid(1:Nqgrid);                % wrap the right boundary

% the second derivative operator
if (matlab>0)
  d2=spdiags([ones(Nqgrid,1) ones(Nqgrid,1) -2*ones(Nqgrid,1) ...
              ones(Nqgrid,1) ones(Nqgrid,1)],[-Nqgrid+1 -1:1 Nqgrid-1], ...
             Nqgrid, Nqgrid)/dx;
else
  global bcsn_overhang
  bcsn_overhang = 25;                   % how long is the wrapping
  d2=zeros(Nqgrid+2*bcsn_overhang,2); 
  d2(:,1) = -2/dx;
  d2(1:end-1,2)=1/dx;
end

% for finding PDFs, not poisson rates, "active" may not exist
if(~exist('bnds')||(length(bnds)==0))
  active = dx*ones(size(qgrid));
else
  active = zeros(size(qgrid));
  %qgrid+0.5dx are edges of bins
  %(L) as well
  [srt ind] = sort([qgrid(:)+0.5*dx; bnds(:,1); bnds(:,2)]);
  srt(srt<0)=0;                         % deal with boundaries overshoots
  srt(srt>L)=L;
  cum=0;
  isactive=0;                           % currently in active state?
  M=size(bnds,1);
  for i=1:length(ind)
    if(ind(i)<=Nqgrid)                  % edge
      if(isactive>0)
        cum = cum + srt(i)-srt(i-1);
      end
      active(ind(i)) = cum;             % add cumulated time
      cum = 0;                          % reset cumulation
    elseif(ind(i)<=Nqgrid+M)            % start active
      isactive = isactive+1;
      if(isactive>1)                    % if the activation state
        cum = cum + srt(i)-srt(i-1);    % was nonzero even before
      end;                              % this time
    else                                % end active
      isactive = isactive-1;
      if(isactive>=0)                   % if was active before this moment
        cum = cum + srt(i)-srt(i-1);
      end
    end
  end
  %last activity wraps around
  active(1)=active(1)+cum;
end

disp(['The fraction of "active" time is: ' num2str(sum(active)/L) ]) 
% now initialize the action
action(0, 0, 0, 0, qgrid, 0, 0, 1);