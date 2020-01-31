function [vm ext_v]=p_DualFloor(v)
% DUALFLOOR determines smallest amount of the combined influence of
% the exercised and preventive power of a pair of coalitions.
%
% Source:  H. I. Meinhardt. The Modiclus Reconsidered. Technical report, Karlsruhe Institute of Technology (KIT), Karlsruhe, Germany,
%          2018. URL http://dx.doi.org/10.13140/RG.2.2.32651.75043.
%
%          Meinhardt (2018), "Analysis of Cooperative Games with Matlab and Mathematica".
%
%          Pelge and Sudhoelter (2007, p.126)
%
% Usage: vm=p_DualFloor(v)
%
% Define variables:
%
%  output:
%  vm        -- Smallest amount of the combined influence of
%               the exercised and preventive power of a pair of coalitions,
%               the so-called dual cover.
%  ext_v.p   -- The primal extension of dual game dv.
%  ext_v.d   -- The dual extension of v.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   03/01/2018        1.0             hme
%


N=length(v);
[~, n]=log2(N);
dv=dual_game(v);


N1=N+1;
n1=2*n;
N2=2^n1-1;
ii=1;
vs1=zeros(N1,N1);
mS=vs1;
parfor k=1:N1
    cl=zeros(N1,1);
    vc=cl;
    for jj=1:N1
          if k>1 && jj >1
             cl(jj)=(k-1)+(jj-1)*N1;
             u=v(k-1)+dv(jj-1);
             w=v(jj-1)+dv(k-1);
             vc(jj)=min(u,w);
           elseif k==1 && jj >1
             cl(jj)=N1*(jj-1);
             u=dv(jj-1);
             w=v(jj-1);
             vc(jj)=min(u,w);
           elseif k>1 && jj==1
             cl(jj)=k-1;
             u=v(k-1);
             w=dv(k-1);
             vc(jj)=min(u,w);
           end
    end
    vs1(k,:)=vc;
    mS(k,:)=cl;
end
clear cl vc v;

vm=zeros(1,N2);
for k=1:N1
   sS=mS(k,:);
   idx=find(sS>0);
   cc=sS(idx);
   vm(cc)=vs1(k,idx);
end
if nargout == 2
  ext_v.p=vs2';
  ext_v.d=vs1';
end
