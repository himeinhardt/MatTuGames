function [cvq A]=p_convex_gameQ(v,tol)
% P_CONVEX_GAMEQ returns 1 whenever the game v is convex using Matlab's PCT.
%
%
% Usage: [cvq A]=p_convex_gameQ(v,tol)
% Define variables:
%  output:
%  cvq      -- Returns 1 (true) or 0 (false).
%  A        -- A pair matrix, which indicates whether the 
%              marginal contribution of the pair (i,j) 
%              has increasing differences (1) or not (0).
%  input:
%  v        -- A TU-game of length 2^n-1.
%  tol      -- Tolerance value. By default, it is set to (-2*10^4*eps).
%              (optional) 


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   05/18/2011        0.1 alpha        hme
%   10/27/2012        0.3              hme
%   05/16/2014        0.5              hme
%                

if nargin<2
   tol=-2*10^4*eps;
end

N=length(v);
[~, n]=log2(N);
A=eye(n);
S=1:N;

parfor i=1:n
  for j=1:n
  if A(i,j)==0
   zi=bitset(0,i); % empty set
   zj=bitset(0,j);
   zij=bitset(zi,j);
   dz=v(zij)-v(zj)-v(zi)>=tol; % empty set 
   Tni=bitget(S,i)==0;
   Tnj=bitget(S,j)==0;
   lg=Tni & Tnj;
   T=S(lg);
   Ti=bitset(T,i);
   Tj=bitset(T,j);
   Tij=bitset(Ti,j);
   dv=v(Tij)-v(Tj)-v(Ti)+v(T)>=tol;
   A(i,j)=all([all(dv),dz]);
  else
  end
 end
end

cvq=all(all(A));
