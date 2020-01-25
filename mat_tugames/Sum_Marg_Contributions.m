function Smarg=Sum_Marg_Contributions(v,T,n,tol)
% Sum_Marg_Contributions returns 1 whenever for coalition T the sum of 
% marginal contributions is positive.
%
% Usage: Smarg=Sum_Marg_Contributions(v,T,n,tol)
% Define variables:
%  output:
%  Smarg    -- Returns 1 whenever the super set T has a positive
%              sum of marginal contributions (1) or a negative (0).
%              Definition of average convexity.
%  input:
%  v        -- A TU-game of length 2^n-1.
%  T        -- A super-set T. It is a positive number.
%  n        -- The number of players involved in game v.
%  tol      -- Tolerance value. By default, it is set to (-2*10^4*eps).
%              (optional) 


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/13/2010        0.1 beta        hme
%                



if nargin<4
   tol=-2*10^4*eps;
end


sT=SubSets(T,n);

lt=length(sT);
tP=cell(lt,1);
Tni=cell(lt,1);
Sni=cell(lt,1);
mrg=cell(lt,1);
zt=zeros(lt,1);
summg=zeros(lt,1);
grSq=false(lt,1);


it=0:-1:1-n;
vecT=rem(floor(sT(:)*pow2(it)),2)==1;


if lt==1
   grSq=v(T)>=tol;
else

 J=1:n;

 for k=1:lt
   tP{k}=J(vecT(k,:));
   Tni{k}=bitset(T,tP{k},0);
   Sni{k}=bitset(sT(k),tP{k},0);
   if Sni{k}==0
     mrg{k}=v(T)-v(Tni{k})-v(sT(k));
     lgm(k)=length(mrg{k});
     summg(k)=mrg{k}*ones(lgm(k),1);
     grSq(k)=summg(k)>=tol;
   else
     mrg{k}=v(T)-v(Tni{k})-v(sT(k))+v(Sni{k});
     lgm(k)=length(mrg{k});
     summg(k)=mrg{k}*ones(lgm(k),1);
     grSq(k)=summg(k)>=tol;
   end
 end
end
Smarg=all(grSq);
