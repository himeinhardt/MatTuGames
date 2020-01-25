function st_x=select_starting_pt(clv)
%SELECT_STARTING_PT computes a starting point for the pre-kernel computation.
% 
%  Usage: st_x=select_starting_pt(clv);
%
%
% Define variables:
%  output:
%  st_x        -- Computed starting point.
%  input:
%  clv         -- TuGame class object.
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   10/29/2012        0.3             hme
%                
    
v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;

S=1:N;
st_x=zeros(1,n);
r_v=zeros(1,n);
Sa=cell(n);

for i=1:n
 a=bitget(S,i)==1;
 Sa{i}=S(a);
 b=Sa{i}-Sa{i}(1);
 bg=b>0;
 s=b(bg);
 vni=[0, v(s)];
 r_v(i)=max(v(Sa{i})-vni);
end

Rm2=r_v(1); for ii=2:n, Rm2=[Rm2 r_v(ii) Rm2+r_v(ii)]; end
Rm=Rm2(N:-1:1);
Rm(1)=[];
Rm(N)=0;

for i=1:n
 vnr=r_v(i)*(v(N)-Rm)./Rm2;
 st_x(i)=max(vnr(Sa{i}));
end

