function [acvQ Smarg]=average_concaveQ(v,tol);
% AVERAGE_CONCAVEQ returns 1 whenever the game v is average.concave.
% For n>14 this function needs some time to complete.
%

% Define variables:
%  output:
%  acvQ     -- Returns 1 (true) or 0 (false).
%  Smarg    -- Returns the list of super-sets which have a negative
%              sum of marginal contributions (1) or a postive (0).
%              The vector has length 2^n-1.
%  input:
%  v        -- A TU-game of length 2^n-1.
%  tol      -- Tolerance value. By default, it is set to 2*10^4*eps.
%              (optional) 


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   06/08/2020       1.9              hme
%                


if nargin<2
   tol=2*10^4*eps;
end


N=length(v);
[~, n]=log2(N);
Smarg=zeros(1,N);
%sumMg=cell(1,N);


for S=1:N;
Smarg(S)=sum_marg_contribution(v,S,n,tol);
end

acvQ=all(Smarg);

%-----------------------------------------
function [Smarg]=sum_marg_contribution(v,T,n,tol);

sT=SubSets(T,n);

lt=length(sT);
tP=cell(lt,1);
Tni=cell(lt,1);
Sni=cell(lt,1);
mrg=cell(lt,1);
summg=zeros(lt,1);
grSq=zeros(lt,1);


subT=dec2bin(sT,n);
flT=fliplr(subT);
vecT=logical(flT-'0');

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
     grSq(k)=summg(k)<=tol;
   else
     mrg{k}=v(T)-v(Tni{k})-v(sT(k))+v(Sni{k});
     lgm(k)=length(mrg{k});
     summg(k)=mrg{k}*ones(lgm(k),1);
     grSq(k)=summg(k)<=tol;
   end
 end
end
Smarg=all(grSq);
