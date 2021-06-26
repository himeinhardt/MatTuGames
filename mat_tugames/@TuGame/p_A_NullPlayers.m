function [anlp,anlpQ]=p_A_NullPlayers(clv,tol)
% P_A_NullPlayers returns the players who are A-Null players using Matlab's PCT.
%
% Usage: adpQ=p_A_NullPlayers(clv,tol)
% Define variables:
%  output:
%  anlp     -- Returns list of A-Null players. Empty set if there is no A-null player.
%  anlpQ    -- Returns 1 (true) whenever the player k is an A-Null player, 
%              otherwise 0 (false). 
%
%  input:
%  clv      -- TuGame class object.
%  tol      -- Tolerance value. Its default value is set to 10^6*eps.
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   02/25/2020        1.9             hme
%                

if nargin<2
 tol=10^6*eps;
end

N=clv.tusize;
n=clv.tuplayers;
v=clv.tuvalues;

S=1:N;
Si=zeros(n,1);
k=1:n;
cli=2.^(k-1);
parfor ss=1:N
   a=k(bitget(ss,k)==1);
   lss=length(a);
   nsi=bitset(ss,k,0);
   nsi=nsi(nsi~=ss);
   if nsi==0
     dvS(ss)=v(ss);
   else	  
     dvS(ss)=sum(v(nsi))/lss;	   
   end
end

anlpQ=false(1,n);
parfor kk=1:n
    cl=S(bitget(S,kk)==1);
    anlpQ(kk)=all(abs(dvS(cl))<tol) && v(cli(kk))==0;
end	
anlp=k(anlpQ);


