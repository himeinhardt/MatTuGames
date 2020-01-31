function shQ=p_ShapleyQ(v,x,tol)
% P_SHAPLEYQ checks if the imputation x is a Shapley value of game v
% using Matlab's PCT.
% 
%  Usage: shQ=p_ShapleyQ(v,x,tol);
%
%
% Define variables:
%  output:
%  shQ      -- Returns 1 (true) whenever the impuatation x is 
%              the Shapley value of v, otherwise 0 (false).
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n) (optional)
%  tol      -- Tolerance value. By default, it is set to 10^6*eps.
%             (optional)


%
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   11/12/2014        0.6             hme
%

if nargin<1
  error('At least a TU game must be specified!');
elseif nargin<2
  x=p_ShapleyValue(v);
  tol=10^6*eps;
elseif nargin<3
  tol=10^6*eps;
end

N=length(v);
[~, n]=log2(N);
if (2^n-1)~=N
    error('Game has not the correct size!');
end

DecG=p_DecomposeGame(v);
shw=p_ShapleyValue(DecG.w);
shz=p_ShapleyValue(DecG.z);
zv=zeros(1,n);
qw=all(abs(zv-shw)<tol);
qz=all(abs(x-shz)<tol);
shQ = qw && qz;
