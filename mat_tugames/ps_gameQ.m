function PS=ps_gameQ(v,tol)
% PS_GAMEQ checks whether the game v is a PS game, if so the
% pre-nucleolus and the Shapley value coincide.
%
% Source: Kar et al. (MSS 57; pp 16-25; 2009) 
%
% Usage: PS=ps_gameQ(clv,tol)
% 
% Define variables:
% output:
%  PS       -- A structure element with the following contents:
%    psQ    -- Returns true (1) or false (0).
%    Q      -- Returns an array of ones and/or zeros of length n.
%    c      -- The sum of the marginal contribution and its
%              complement for each player. 
%    Mi     -- The matrix of marginal contribution and its
%              complement for each player. Each row must be constant. 
%    sh     -- The Shapley value of the PS game.
%    pnc    -- The pre-nucleolus which must be equal to the Shapley value.
%              It is (1/2)*c. 
%
%
% input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  tol      -- Tolerance value. Its default value is set to 10^6*eps.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   01/25/2015        0.6             hme
%

if nargin < 2
   tol=10^6*eps;
end

N=length(v);
[~, n]=log2(N);
S=1:N;
k=1:n;
iC=2.^(k-1);
li=2^(n-1);
Mi=zeros(n,li);
dv=dual_game(v);
c=zeros(1,n);

for ii=1:n 
    Si = S(bitget(S,ii)==1);
    Swi = Si-iC(ii);
    mi0=v(Si(1))-0;
    dmi0=dv(Si(1))-0;
    Si(1)=[];
    Swi(1)=[];
    mi=v(Si)-v(Swi);
    dmi=dv(Si)-dv(Swi);
    Mi(ii,:)=[mi0+dmi0,mi+dmi];
    mc=max(Mi(ii,:));
    psQ(ii)=all(abs(Mi(ii,:)-mc)<tol);
    if psQ(ii)==1
       c(ii)=mc;
    else
       c(ii)=-inf; 
    end    
end

PS.psQ=all(psQ);
PS.Q=psQ;
PS.c=c;
PS.Mi=Mi;
if PS.psQ==1
   PS.sh=(1/2)*c;
   PS.pnc=PS.sh;
else
   PS.sh=ShapleyValue(v);
   try
      PS.pnc=cplex_prenucl_llp(v);
   catch
      PS.pnc=PreNucl(v);
   end
end    
