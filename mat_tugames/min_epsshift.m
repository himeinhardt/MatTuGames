function V_MEPS=min_epsshift(v)
% MIN_EPSSHIFT computes for an almost-convex game the min epsilon shift to construct a convex game. 
% More general approach as by  J. Getan et. al.
%
% Source: J. Getan et. al. (2012), The bargaining set and the kernel for almost-convex games. 
%
% Define structure variables:
%  output:
%  v_meps   -- Constructed strong epsilon game from value meps.
%  meps     -- Returns the possible min epsilon shift to get a convex game.
%  cv_mepsQ -- Returns 1 (true) whenever game v_meps is convex, otherwise 0. 
%  alcvQ    -- Returns 1 (true) whenever game v is almost-convex, otherwise 0.
%  A        -- A pair matrix, which provides the difference of 
%              marginal contributions of the pair (i,j).
%  v        -- original TU game (input).
%  input:
%  v        -- A TU-game of length 2^n-1.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   11/27/2020        1.9             hme
%                

if nargin<1
   error('At least a TU game must be given as input argument!');
else
   N=length(v);
   [~, n]=log2(N);	
   if (2^n-1)~=N
      error('Game has not the correct size!');
   end
end

A=diag(-inf(1,n));
S=1:N;
Tni=cell(n);
Tnj=cell(n);
Ti=cell(n);
Tj=cell(n);
Tij=cell(n);
T=cell(n);
lg=cell(n);
dv=cell(n);
zij=eye(n);
dz=eye(n);

for i=1:n
  for j=1:n
  if A(i,j)==0
   zi(i,j)=bitset(0,i); % empty set
   zj(i,j)=bitset(0,j);
   zij(i,j)=bitset(zi(i,j),j);
   dz(i,j)=v(zj(i,j))+v(zi(i,j))+v(zij(i,j)); % empty set 
   Tni{i,j}=bitget(S,i)==0;
   Tnj{i,j}=bitget(S,j)==0;
   lg{i,j}=Tni{i,j} & Tnj{i,j};
   T{i,j}=S(lg{i,j});
   Ti{i,j}=bitset(T{i,j},i);
   Tj{i,j}=bitset(T{i,j},j);
   Tij{i,j}=bitset(Ti{i,j},j);
   dv{i,j}=v(Tj{i,j})+v(Ti{i,j})-(Tij{i,j})-v(T{i,j});
   A(i,j)=min([min(dv{i,j}),dz(i,j)]);
  else
  end
 end
end
 
meps=max(max(A));
v_meps=streps_value(v,meps);
alcvQ=AlmostConvex_gameQ(v);
cv_mepsQ=convex_gameQ(v_meps);
V_MEPS.meps=meps;
V_MEPS.cv_mepsQ=cv_mepsQ;
V_MEPS.alcvQ=alcvQ;
V_MEPS.v_meps=v_meps;
V_MEPS.A=A;
V_MEPS.v=v;
