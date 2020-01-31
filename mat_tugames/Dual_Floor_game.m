function v_t=Dual_Floor_game(v,x)
% DUAL_Floor_GAME computes from (v,x) a modified Davis-Maschler anti-reduced game vS on S at x for
% game v.
%
% Usage: vt=Dual_Floor_game(v,x)
% Define variables:
%  output:
%  v_t     -- A set of modified Davis-Maschler anti-reduced game vS w.r.t. x.
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
%   02/06/2018        0.9             hme
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
v_t=cell(1,N-1);
%dv_t=cell(2,N-1);


vS=DM_AntiReduced_game(v,x);
dvS=DM_AntiReduced_game(dv,x);


for k=1:N-1
    v_t{1,k}=min(vS{1,k}+mx_v+2*mx_dv,dvS{1,k}+mx_dv+2*mx_v);
    v_t{1,k}(end)=vS{1,k}(end);
end
