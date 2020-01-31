function [v_sp, mat_hd, MatW, MatV, A, mat_vz]=game_space(clv,x,slc,smc)
% GAME_SPACE computes the game space which replicates x as a pre-kernel element.
%
%
% Usage: [v_sp mat_hd MatW MatV A mat_vz]=clv.game_space(x,slc,smc)
%
% Define variables:
% output:
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
%   05/10/2014        0.5              hme
%                

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;

onm=true(n);
upe=tril(onm,-1);
[~, A]=clv.BestCoalitions(x,smc);
trA=A';
drij=trA(upe)';
drji=A(upe)';
uG=speye(N); % unity games

MatV=uG(drji,:)-uG(drij,:);
MatV(end+1,:)=uG(N,:);
MatV=sparse(MatV);
[hd, gb]=clv.unanimity_games();
MatW=MatV*gb;
MatW=full(MatW);
nlW=null(MatW);
sW=size(nlW);
hd=hd';
HDm=repmat(hd,1,sW(2));
mat_hz=slc*nlW;
mat_vz=gb*mat_hz;
mat_vz=mat_vz';
mat_hd=HDm+mat_hz;
w_sp=gb*mat_hd;
v_sp=w_sp';


