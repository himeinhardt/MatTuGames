function v_t=Dual_Cover_game(clv,x)
% DUAL_COVER_GAME computes from (v,x) a modified Davis-Maschler reduced game vS on S at x for
% game v.
%
% Usage: vt=clv.Dual_Cover_game(x)
% Define variables:
%  output:
%  v_t     -- A set of modified Davis-Maschler reduced game vS w.r.t. x.
%  input:
%  clv      -- TuGame class object.
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


v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
dv=clv.dual_game();

exc_v=clv.excess(x);
exc_dv=excess(dv,x);
[sx_v,idx_v]=sort(exc_v,'descend');
[sx_dv,idx_dv]=sort(exc_dv,'descend');
mx_v=sx_v(1);
mx_dv=sx_dv(1);
v_t=cell(1,N-1);
%dv_t=cell(2,N-1);


vS=clv.DM_Reduced_game(x);
dvS=DM_Reduced_game(dv,x);


for k=1:N-1
    v_t{1,k}=max(vS{1,k}+mx_v+2*mx_dv,dvS{1,k}+mx_dv+2*mx_v);
    v_t{1,k}(end)=vS{1,k}(end);
end
