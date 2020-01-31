function v_x=p_ECFloorGame(v,x)
% P_ECFloorGame computes from (v,x) an excess comparability floor of game v using Matlab's PCT.
%
% Source:  H. I. Meinhardt. The Modiclus Reconsidered. Technical report, Karlsruhe Institute of Technology (KIT), Karlsruhe, Germany,
%          2018. URL http://dx.doi.org/10.13140/RG.2.2.32651.75043.
%
%          Meinhardt (2018), "Analysis of Cooperative Games with Matlab and Mathematica".
%
%
% Usage: vt=p_ECFloorGame(v,x)
%
% Define variables:
%  output:
%  v_x     -- The excess comparability floor game w.r.t. x.
%
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
sx_v=sort(exc_v,'ascend');
sx_dv=sort(exc_dv,'ascend');
mx_v=sx_v(1);
mx_dv=sx_dv(1);
v_x=zeros(1,N);




parfor k=1:N-1
    v_x(k)=min(v(k)+mx_v+2*mx_dv,dv(k)+mx_dv+2*mx_v);
end
v_x(N)=v(N);
