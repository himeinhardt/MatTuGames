function v_t=p_DM_Anti_Derived_game(v,x)
% p_DM_ANTI_DERIVED_GAME computes from (v,x) a modified Davis-Maschler anti-reduced game vS on S at x for
% game v using Matlab's PCT.
%
% Usage: vt=p_DM_Anti_Derived_game(v,x)
% Define variables:
%  output:
%  v_t      -- A set of modified Davis-Maschler anti-reduced game vS w.r.t. x.
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n).
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

n=length(x);
N=length(v);

exc_v=excess(v,x);
dv=dual_game(v);
exc_dv=excess(dv,x);
[sx_v,idx_v]=sort(exc_v,'ascend');
[sx_dv,idx_dv]=sort(exc_dv,'ascend');
mx_v=sx_v(1);
mx_dv=sx_dv(1);
v_t=cell(1,N);
%dv_t=cell(2,N-1);


vSa=p_DM_AntiReduced_game(v,x);
vS=vSa{:,1};
clear vSa;
dvSa=p_DM_AntiReduced_game(dv,x);
dvS=dvSa{:,1};
clear dvSa;

parfor k=1:N-1
    mv=min(vS{1,k}-mx_v,dvS{1,k}-mx_dv);
    mv(end)=vS{1,k}(end);
    v_t{1,k}=mv;
end
v_t{1,N}=min(v-mx_v,dv-mx_dv);
v_t{1,N}(end)=v(N);
