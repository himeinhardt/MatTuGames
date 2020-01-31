function [v_sp, A]=p_game_space_red(clv,x,slc,smc)
% P_GAME_SPACE_RED computes the game space which replicates x as a pre-kernel element.
% Same as the function p_game_space() but with less output arguments to save memory.
% Using Matlab's PCT.
%
% Usage: [v_sp A]=clv.p_game_space(x,slc,smc)
%
% Define variables:
% output:
%  v_spc             -- Game space spanned by the basis of the null space
%                       MatW.
%  A                 -- Indicates the set of equivalence classes/most
%                       effective coalitions w.r.t. the pre-kernel
%                       element x.
%
%
%  input:
%  clv    -- TuGame class object.
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
%   10/29/2012        0.3              hme
%   05/15/2014        0.5              hme
%                

N=clv.tusize;
n=clv.tuplayers;

onm=true(n);
upe=tril(onm,-1);
[~, A]=clv.p_BestCoalitions(x,smc);
trA=A';
drij=trA(upe)';
drji=A(upe)';
uG=speye(N); % unity games

MatV=uG(drji,:)-uG(drij,:);
MatV(end+1,:)=uG(N,:);
MatV=sparse(MatV);
[hd, gb]=clv.p_unanimity_games();
MatW=MatV*gb;
clear MatV;
MatW=full(MatW);
nlW=null(MatW);
clear MatW;
[~,sW2]=size(nlW);
hd=hd';
HDm=repmat(hd,1,sW2);
mat_hz=slc*nlW;
clear nlW;
mat_hd=HDm+mat_hz;
clear mat_hz;
if n>13
 v_sp=zeros(sW2,N);
 parfor k=1:sW2
   v_sp(k,:)=gb*mat_hd(:,k);
 end
 clear mat_hd gb;
else
 spmd
   cgb=codistributed(double(gb));
   cmhd=codistributed(mat_hd);
   v_sp1=cgb*cmhd;
 end
 clear mat_hd gb;
 v_sp=gather(v_sp1)';
end
