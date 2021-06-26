function APDP=AP_DummyPlayer_propertyQ(clv,x,tol)
% AP_DummyPlayer_propertyQ checks if the solution x satisfies the AP-Dummy player property.
%
% Source: Emilio Calvo Ramón and Esther Gutiérrez-López (2020), 
%         The Equal Collective Gains value in Cooperative Games
%
%
% Usage: APDP=clv.AP_DummyPlayer_propertyQ(x,tol)
% Define variables:
%  output:
%  Q        -- Returns 1 (true) whenever the solution satisfies the AP-Dummy player property, 
%              otherwise 0 (false). 
%  ply      -- Returns 1 (true) whenever the player k is an AP-Dummy player, 
%              otherwise 0 (false).
%  x        -- Returns input payoff vector.
%
%  input:
%  clv      -- TuGame class object.
%  x        -- payoff vector of size(1,n). Must be efficient.
%  tol      -- Tolerance value. Its default value is set to 10^6*eps.
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   12/25/2020        1.9             hme
%                

if nargin<3
 tol=10^7*eps;
end

N=clv.tusize;
n=clv.tuplayers;
v=clv.tuvalues;

effQ=abs(sum(x)-v(N))<tol;

if effQ==0
   APDP.Q=false;
   APDP.ply=[];
   APDP.x=x;
   return;
end


S=1:N;
Si=zeros(n,1);
k=1:n;
cli=2.^(k-1);
for ss=1:N
   a=k(bitget(ss,k)==1);
   lss=length(a);
   nsi=bitset(ss,k,0);
   nsi=nsi(nsi~=ss);
   if nsi==0
     dvi(ss)=v(ss);
   else	   
     dvi(ss)=lss*v(ss)-sum(v(nsi));
   end
   dvS(ss)=(dvi(ss)-v(ss))/lss;
end

APpQ=false(1,n);
apdpQ=false(1,n);
for kk=1:n
    cl=S(bitget(S,kk)==1);
    apdpQ(kk)=all(abs(dvS(cl))<tol);
    if apdpQ(kk)==1
       APpQ(kk)=x(kk)==v(cl(kk));
    elseif x(kk)~=v(cl(kk)) %% contrapositive
       APpQ(kk)=apdpQ(kk)==0;	    
    end	    
end	
APDP.Q=all(APpQ);
APDP.ply=apdpQ;
APDP.x=x;

