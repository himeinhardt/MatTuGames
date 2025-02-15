function eS=essentialSet(clv,tol)
% ESSENTIALSET computes the set of essential coalitions.
%
%
% It requires Partition of an Integer from
% 
% SOURCE: http://www.mathworks.com/matlabcentral/fileexchange/12009-partitions-of-an-integer
%
% Usage: eS=clv.essentialSet()
%
% 
% Define variables:
% output:
% eS        -- Collection of essential coalitions.
%
% input: 
%  clv        -- TuGame class object.
%  tol        -- Tolerance value. Its default value is set to 10^6*eps.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/11/2017        0.9             hme
%                

if nargin<2
 tol=10^6*eps;
end

v=clv.tuvalues;
n=clv.tuplayers;

jj=1:n;
eS=[];
lrQ=[];
it=0:-1:1-n;
si=2.^(jj-1);
for S=1:N
   pt=intpartitions(S);
   sm=numel(pt);
   if any(S==si)==0
      for k=1:sm-1
          ui=unique(pt{k})';
          lui=length(ui);
          sQ=sum(ui);
          if sQ==S
             ps=bitget(S,jj)==1;
             lp=ps*ones(n,1);
             sS=double(ui);
             ci=rem(floor(sS(:)*pow2(it)),2)==1;
             sci=size(ci,1);
             ad=ones(1,sci)*ci;
             if all(ad==ps) 
                vSQ=v(S)>sum(v(ui))+tol;
                lrQ=[vSQ,lrQ];
             end     
          end
      end
   else
     lrQ=[true,lrQ]; 
   end    
   if all(lrQ)
      eS=[eS,S];
   end
   lrQ=[];
end
