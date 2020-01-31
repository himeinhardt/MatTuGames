function [v_sp, A]=p_game_space_red(v,x,slc,smc,method)
% P_GAME_SPACE_RED computes the game space which replicates x as a pre-kernel element.
% Same as the function p_game_space() but with less output arguments to save memory.
% Using Matlab's PCT.
%
% Usage: [v_sp A]=p_game_space(v,x,slc,smc)
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
%  v                 -- A Tu-Game v of length 2^n-1. 
%  x                 -- pre-kernel payoff of length(1,n)
%  scl               -- scaling factor
%  smc               -- selecting from effc the smallest/largest 
%                       cardinality (optional). Value 1 or 0.
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
%   06/02/2011        0.1 alpha        hme
%   05/15/2014        0.5              hme
%                

if nargin < 5
   method = 'full';
end


N=length(v); n=length(x);
onm=true(n);
upe=tril(onm,-1);
[~, A]=p_BestCoalitions(v,x,smc);
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
[hd, gb]=p_unanimity_games(v);
MatW=MatV*gb;
clear MatV uG;
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
