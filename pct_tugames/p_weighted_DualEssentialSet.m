function wdeS=p_weighted_DualEssentialSet(v,pS,tol)
% P_WEIGHTED_DUALESSENTIALSET computes the set of weighted dually essential
% coalitions using Matlab's PCT. For an alternative computational approach see the
% function weighted_DuallyEssentialSet().
%
% It requires Partition of an Integer from
% 
% SOURCE: http://www.mathworks.com/matlabcentral/fileexchange/12009-partitions-of-an-integer
%
% Usage: wdeS=weighted_DualEssentialSet(v,pS)
%
% 
% Define variables:
% output:
% eS        -- Collection of weighted dually essential coalitions.
%
% input: 
%  v        -- A Tu-Game v of length 2^n-1. 
%  pS       -- A vector of weights of length 2^n-1. Default is per capita weights.
%  tol      -- Tolerance value. Its default value is set to 10^6*eps.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/10/2017        0.9             hme
%                


if nargin<2
 pS='';   
 tol=10^6*eps; % Change this value if the solution is not correct.
elseif nargin<3
 tol=10^6*eps;   
end


N=length(v);
[~, n]=log2(N);

if isempty(pS)
   sS=1:N;
   for k=1:n, A1(:,k) = -bitget(sS,k);end
   mat=-A1';
   clS=ones(1,n)*mat;
   pS=1./clS;
   pS(N)=1;
   clear sS;
end

pv=pS.*v;

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
   if any(NS==si)==0 && S>0
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
                ti=N-ui;
                vSQ=pv(S)>sum(pv(ti))-(lui-1)*pv(N)+tol;
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
      wdeS=[wdeS,S];
   end
%   lrQ=[];
end
cnt
wdeS(end+1)=N;

