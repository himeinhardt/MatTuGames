function [cvq A]=convex_gameQ(v,tol)
% CONVEX_GAMEQ returns 1 whenever the game v is convex.
%
%
% Usage: [cvq A]=convex_gameQ(v,tol)
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
%   08/06/2010        0.1 beta        hme
%   10/27/2012        0.3             hme
%                

if nargin<2
   tol=-2*10^4*eps;
end

N=length(v);
[~, n]=log2(N);
A=eye(n);
S=1:N;
Tni=cell(n);
Tnj=cell(n);
Ti=cell(n);
Tj=cell(n);
Tij=cell(n);
T=cell(n);
lg=cell(n);
dv=cell(n);
zi=eye(n);
zj=eye(n);
zij=eye(n);
dz=eye(n);

for i=1:n
  for j=1:n
  if A(i,j)==0
   zi(i,j)=bitset(0,i); % empty set
   zj(i,j)=bitset(0,j);
   zij(i,j)=bitset(zi(i,j),j);
   dz(i,j)=v(zij(i,j))-v(zj(i,j))-v(zi(i,j))>=tol; % empty set 
   Tni{i,j}=bitget(S,i)==0;
   Tnj{i,j}=bitget(S,j)==0;
   lg{i,j}=Tni{i,j} & Tnj{i,j};
   T{i,j}=S(lg{i,j});
   Ti{i,j}=bitset(T{i,j},i);
   Tj{i,j}=bitset(T{i,j},j);
   Tij{i,j}=bitset(Ti{i,j},j);
   dv{i,j}=v(Tij{i,j})-v(Tj{i,j})-v(Ti{i,j})+v(T{i,j})>=tol;
   A(i,j)=all([all(dv{i,j}),dz(i,j)]);
  else
  end
 end
end
 
cvq=all(all(A));
