function st_x=p_select_starting_pt(v)
%P_SELECT_STARTING_PT computes a starting point for the pre-kernel computation
% while using Matlab's PCT.
% 
%  Usage: st_x=p_select_starting_pt(v);
%
%
% Define variables:
%  output:
%  st_x        -- Computed starting point.
%  input:
%  v           -- A Tu-Game v of length 2^n-1. 
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   06/04/2012        0.2 beta        hme
%   09/14/2012        0.2             hme
%   10/27/2012        0.3             hme
%                
    
if nargin<1
    error('The game must be given!');
else
   N=length(v);
   [~, n]=log2(N);
end

st_x=zeros(1,n);
r_v=zeros(1,n);
s=cell(n,1);
Sa=cell(n,1);
S=1:N;

parfor i=1:n
 a=bitget(S,i)==1;
 Sa=find(a);
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
mr=max((v(N)-Rm)./Rm2);

k=1:n;
st_x(k)=r_v(k)*mr;

