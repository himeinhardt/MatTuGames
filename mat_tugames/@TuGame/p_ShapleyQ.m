function shQ=p_ShapleyQ(clv,x,tol)
% P SHAPLEYQ checks if the imputation x is a Shapley value of game v
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
%  clv      -- TuGame class object.
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

if nargin<2
  x=clv.p_ShapleyValue();
  tol=10^6*eps;
elseif nargin<3
  tol=10^6*eps;
end

N=clv.tusize;
n=clv.tuplayers;

DecG=clv.p_DecomposeGame;
shw=p_ShapleyValue(DecG.w);
shz=p_ShapleyValue(DecG.z);
zv=zeros(1,n);
qw=all(abs(zv-shw)<tol);
qz=all(abs(x-shz)<tol);
shQ = qw && qz;
