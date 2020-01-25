function [v_sp A]=p_game_space_red(v,x,slc,smc)
% P_GAME_SPACE_RED computes the game space which replicates x as a pre-kernel element.
% Same as the function p_game_space() but with less output arguments to save memory.
% Using Matlab's PCT.
%
% Usage: [v_sp A]=p_game_space(v,x,slc,smc)
% Define variables:
% output:
%  v_spc             -- Game space spanned by the basis of the null space
%                       MatW.
%  A                 -- Indicates the set of equivalence class/most
%                       effective coalitions w.r.t. the pre-kernel
%                       element x.
%
%
%  input:
%  v      -- A Tu-Game v of length 2^n-1. 
%  x      -- pre-kernel payoff of length(1,n)
%  scl    -- scaling factor
%  smc    -- selecting from effc the smallest/largest 
%            cardinality (optional). Value 1 or 0.

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   06/02/2011        0.1 alpha        hme
%                



N=length(v); n=length(x);
S=1:N;
onm=ones(n);
drij=zeros(1,n);
drji=zeros(1,n);
upe=logical(tril(onm,-1));
[e A]=p_BestCoalitions(v,x,smc);
trA=A';
drij=trA(upe)';
drji=A(upe)';
uG=eye(N); % unity games

MatV=uG(:,drji)-uG(:,drij);
MatV(:,end+1)=uG(:,N);
MatV=sparse(MatV);
alpvec=MatV'*v';
[uc gb]=p_unanimity_games(v);
MatW=MatV'*gb;
clear MatV;
MatW=full(MatW);
nlW=null(MatW);
clear MatW;
sW=size(nlW);
hd=uc';
HDm=repmat(hd,1,sW(2));
mat_hz=slc*nlW;
clear nlW;
mat_hd=HDm+mat_hz;
clear mat_hz;
w_sp=gb*mat_hd;
clear mat_hd;
v_sp=w_sp';
