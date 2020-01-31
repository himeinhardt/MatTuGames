function [v_sp, mat_hd, MatW, MatV, A, mat_vz]=game_space(v,x,slc,smc,method)
% GAME_SPACE computes the game space which replicates x as a pre-kernel element.
%
%
% Usage: [v_sp mat_hd MatW MatV A mat_vz]=game_space(v,x,slc,smc)
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
%  input:
%  v                 -- A Tu-Game v of length 2^n-1. 
%  x                 -- pre-kernel payoff of length(1,n)
%  scl               -- scaling factor
%  smc               -- selecting from effc the smallest/largest 
%                        cardinality (optional). Value 1 or 0.
%  method            -- Permissible methods are:
%                       'sparse'  or 
%                       'full'
%                       to compute spares or dense matrices. Default is full.
%                       If you have a memory problem use 'sparse' instead.

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   02/02/2011        0.1 beta        hme
%   05/10/2014        0.5             hme
%                

if nargin < 5
   method = 'full';
end

N=length(v); n=length(x);
onm=true(n);
upe=tril(onm,-1);
[~, A,~]=BestCoalitions(v,x,smc);
trA=A';
drij=trA(upe)';
drji=A(upe)';
if strcmp(method,'full')
  uG=eye(N);
else
  uG=speye(N); % unity games
end
MatV=uG(drji,:)-uG(drij,:);
MatV(end+1,:)=uG(N,:);
[hd, gb]=unanimity_games(v);
MatW=MatV*gb;
if issparse(MatW)
   try 
%  Returns a sparse orthonormal basis for the null space
%  from the QR decomposition.
     nlW=spnull(MatW);
   catch
% SVD decomposition.
     MatW=full(MatW);
     nlW=null(MatW);
   end
else
   nlW=null(MatW);
end
sW=size(nlW);
hd=hd';
HDm=repmat(hd,1,sW(2));
mat_hz=slc*nlW;
mat_vz=(gb*mat_hz)';
mat_hd=HDm+mat_hz;
v_sp=(gb*mat_hd)';

