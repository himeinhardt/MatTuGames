function DPPQ=DummyPlayer_propertyQ(v,x,tol)
% DUMMYPLAYER_PROPERTYQ checks if the solution x satisfies the dummy player property. 
%
%  Usage: DPPQ=DummyPlayer_propertyQ(v,x,tol)
%
% Define variables:
%  output: Fields
%  DppQ     -- Returns true (1), if x satisfies the dummy player property, otherwise
%              false (0).
%  sdp      -- Set of dummy players.
%  pzp      -- Players with zero payoff.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n). Must be efficient.
%  tol      -- Tolerance value. By default, it is set to 10^6*eps.
%              (optional) 
%              

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   12/27/2020        1.9             hme
%

if nargin<2
 x=ShapleyValue(v);   
 tol=10^6*eps;
elseif nargin<3 
 tol=10^6*eps;    
end

N=length(v);
[~, n]=log2(N);
pl=1:n;
cli=2.^(pl-1);
S=1:N;
npQ=false(1,n);
NppQ=false;
for k=1:n
    a=bitget(S,k)==1;
    Swk=S(a==0);
    Sk=S(a);
    vS=[0,v(Swk)];
    npQ(k)=all(abs(v(Sk)-vS-v(cli(k)))<tol);
end    
NpQ=any(npQ);
si=cli(npQ);
snp=pl(npQ);
nnp=pl(npQ==0);
nlp=pl(x(pl)~=v(cli));
if NpQ==1
   zpQ=abs(x(snp)-v(si))<tol;
   pzp=pl(zpQ);
   NppQ=any(zpQ);
elseif any(nlp) %% contrapositive
   slc=ismember(nnp,nlp);	
   NppQ=all(nnp(slc)==nlp);
   pzp=pl(x(pl)==v(cli));
else
   NppQ=false;
   pzp=[];
end

DPPQ=struct('DppQ',NppQ,'sdp',snp,'pzp',pzp);

