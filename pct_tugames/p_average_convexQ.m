function [acvQ Smarg sumMg]=p_average_convexQ(v,tol)
% P_AVERAGE_CONVEXQ returns 1 whenever the game v is average-convex
% using Matlab's PCT.
%
% Usage: [acvQ Smarg sumMg]=p_average_convexQ(v,tol)
% Define variables:
%  output:
%  acvQ     -- Returns 1 (true) or 0 (false).
%  Smarg    -- Returns the list of super-sets which have a positive
%              sum of marginal contributions (1) or a negative (0).
%              The vector has length 2^n-1.
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
%   05/19/2011        0.1 alpha        hme
%   06/27/2012        0.2 beta         hme
%   10/27/2012        0.3              hme
%   05/16/2014        0.5              hme
%   08/02/2018        1,0             hme
%                


if nargin<2
   tol=-2*10^4*eps;
end


N=length(v);
[~, n]=log2(N);
Smarg=zeros(1,N);
sumMg=cell(1,N);

parfor S=1:N;
[Smarg(S) sumMg{S}]=sum_marg_contribution(v,S,n,tol);
end

acvQ=all(Smarg);

%-----------------------------------------
function [Smarg summg]=sum_marg_contribution(v,T,n,tol);

sT=SubSets(T,n);

lt=length(sT);
tP=cell(lt,1);
Tni=cell(lt,1);
Sni=cell(lt,1);
mrg=cell(lt,1);
summg=zeros(lt,1);
grSq=false(lt,1);


it=0:-1:1-n;
vecT=rem(floor(sT(:)*pow2(it)),2)==1;

if lt==1 %% S is empty or S=T
   %grSq=v(T)>=tol;
   grSq=true;
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
