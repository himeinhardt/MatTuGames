function v_t=HMS_Derived_game(v,x,str)
% DM_DERIVED_GAME computes from (v,x,S) a modified Hart-Mas-Colell reduced game vS on S at x for
% game v.
%
% Source:  H. I. Meinhardt. The Modiclus Reconsidered. Technical report, Karlsruhe Institute of Technology (KIT), Karlsruhe, Germany,
%          2018. URL http://dx.doi.org/10.13140/RG.2.2.32651.75043.
%
%          Meinhardt (2018), "Analysis of Cooperative Games with Matlab and Mathematica".
%
%
%
% Usage: vt=HMS_Derived_game(v,x)
%
% Define variables:
%
%  output:
%  v_t     -- A set of modified Davis-Maschler reduced game vS w.r.t. x.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n).
%  str      -- A string that defines different Methods. 
%              Permissible methods are: 
%              'PRN' that is, the Hart-MasColell derived game 
%               in accordance with the pre-nucleolus.
%              'MODIC' that is, the Hart-MasColell derived game 
%               equivalence in accordance with the modiclus.
%              'PRK' that is, the Hart-MasColell derived game 
%               in accordance with pre-kernel solution.
%              'SHAP' that is, Hart-MasColell derived game 
%               in accordance with the Shapley Value, 
%               which is, the original definition. 
%              Default is 'SHAP'.
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
[sx_v,idx_v]=sort(exc_v,'descend');
[sx_dv,idx_dv]=sort(exc_dv,'descend');
mx_v=sx_v(1);
mx_dv=sx_dv(1);
v_t=cell(1,N-1);
%dv_t=cell(2,N-1);


vS=HMS_Reduced_game(v,x,str);
dvS=HMS_Reduced_game(dv,x,str);


for k=1:N-1
    v_t{1,k}=max(vS{1,k}-mx_v,dvS{1,k}-mx_dv);
    v_t{1,k}(end)=vS{1,k}(end);
end
v_t{1,N}=max(v-mx_v,dv-mx_dv);
v_t{1,N}(end)=v(N);
