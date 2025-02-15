function deS=p_DuallyEssentialSet(v,tol)
% DUALLYESSENTIALSET computes the set of dually essential
% coalitions using Matlab's PCT. For an alternative computational approach see the
% function DualEssentialSet().
%
% It requires Partition of an Integer from
% 
% SOURCE: http://www.mathworks.com/matlabcentral/fileexchange/12009-partitions-of-an-integer
%
% Usage: deS=p_DuallyEssentialSet(v)
%
% 
% Define variables:
% output:
% deS        -- Collection of dually essential coalitions.
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
%   08/09/2017        0.9             hme
%                

if nargin<2
 tol=10^6*eps;
end


N=length(v);
[~, n]=log2(N);
dv=dual_game(v);
jj=1:n;
deS=[];
%lrQ=[];
it=0:-1:1-n;
si=2.^(jj-1);
Ni=N-si;
parfor S=1:N-1
   lrQ=[];
   NS=N-S;
   pt=intpartitions(NS);
   sm=numel(pt);
   if any(NS==si)==0 
      for k=1:sm-1
          ui=unique(pt{k})';
          lui=length(ui);
          sQ=sum(ui);
          if sQ==NS
             ps=bitget(NS,jj)==1;
             lp=ps*ones(n,1);
             sS=double(ui);
             ci=rem(floor(sS(:)*pow2(it)),2)==1;
             sci=size(ci,1);
             ad=ones(1,sci)*ci;
             if all(ad==ps)
                vSQ=dv(NS)<sum(dv(ui))-tol;
                lrQ=[vSQ,lrQ];
                if vSQ==0
                   continue;
                end
             end     
          end
      end
   else
       lrQ=[true,lrQ]; 
   end    
   if all(lrQ)
      deS=[deS,S];
   end
%   lrQ=[];
end
deS(end+1)=N;
