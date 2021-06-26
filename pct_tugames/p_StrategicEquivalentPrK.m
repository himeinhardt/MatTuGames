function [x1,y]=p_StrategicEquivalentPrK(v,x)
% P_STRATEGICEQUIVALENTPRK computes a pre-kernel element of game v while referring on 
% a strategic equivalent game using Matlab's PCT.
%
%  Usage: [x1,y]=p_StrategicEquivalentPrK(v,x)
%
% Define variables:
%  output:
%  x1       -- Pre-Kernel element (output) of game v.
%  y        -- Pre-Kernel of a zero-normalized game zv. 
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n) (optional)
%

%
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   10/31/2015        0.7             hme
%

N=length(v);
[~, n]=log2(N);
k=1:n;
ci=2.^(k-1);


zv=zero_normalization(v);

if nargin < 2
   y=p_PreKernel(zv);
else
   y=p_PreKernel(zv,x);
end
x1=y+v(ci);


