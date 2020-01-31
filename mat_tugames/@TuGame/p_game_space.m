function [v_sp, mat_hd, MatW, MatV, A, mat_vz]=p_game_space(clv,x,slc,smc)
% P_GAME_SPACE computes the game space which replicates x as a pre-kernel element.
% Using Matlab's PCT.
%
% Usage: [v_sp mat_hd MatW MatV A mat_vz]=clv.p_game_space(x,slc,smc)
%
% Define variables:
% output:
%  v_sp              -- Game space spanned by the basis of the null space
%                       MatW.
%  mat_hd            -- basis of the null space MatW.
%  MatW              -- Matrix given by equation 7.16 on page 152 Meinhardt
%                       2013.
%  MatV              -- Matrix obtained by the unity games that
%                       constitutes balancedness of the maximal surpluses.
%  A                 -- Matrix of most effective coalitions (eq. class).
%  mat_vz            -- Set of null games.  
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

v=clv.tuvalues;
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
MatW=full(MatW);
nlW=null(MatW);
[~,sW2]=size(nlW);
hd=hd';
HDm=repmat(hd,1,sW2);
mat_hz=slc*nlW;
clear nlW hd;
if n>13
 mat_vz=zeros(sW2,N);
 parfor k=1:sW2
   mat_vz(k,:)=gb*mat_hz(:,k);
 end
else
 spmd
   cmhz=codistributed(mat_hz);
   cgb=codistributed(double(gb));
   mat_vz1=cgb*cmhz;
 end
 mat_vz=gather(mat_vz1)';
 clear mat_vz1;
end
mat_hd=HDm+mat_hz;
clear HDm mat_hz mat_vz1;
if n>13
 v_sp=zeros(sW2,N);
 parfor k=1:sW2
   v_sp(k,:)=gb*mat_hd(:,k);
 end
 clear gb;
else
 spmd
   cmhd=codistributed(mat_hd);
   v_sp1=cgb*cmhd;
 end
 v_sp=gather(v_sp1)';
end
