function st_x=p_select_starting_pt(clv)
%P_SELECT_STARTING_PT computes a starting point for the pre-kernel computation
% while using Matlab's PCT.
% 
%  Usage: st_x=p_select_starting_pt(clv);
%
%
% Define variables:
%  output:
%  st_x        -- Computed starting point.
%
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
%   10/30/2012        0.3             hme
%                
    
v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
S=1:N;
st_x=zeros(1,n);
r_v=zeros(1,n);

parfor i=1:n
 a=bitget(S,i)==1;
 Sa=S(a);
 b=Sa-Sa(1);
 bg=b>0;
 s=b(bg);
 vni=[0, v(s)];
 r_v(i)=max(v(Sa)-vni);
end

Rm2=r_v(1); for ii=2:n, Rm2=[Rm2 r_v(ii) Rm2+r_v(ii)]; end
Rm=Rm2(N:-1:1);
Rm(1)=[];
Rm(N)=0;

parfor i=1:n
 a=bitget(S,i)==1;
 vnr=r_v(i)*(v(N)-Rm)./Rm2;
 st_x(i)=max(vnr(a));
end

