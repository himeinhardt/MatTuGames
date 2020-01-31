function [RepSol, RepBMat]=p_replicate_prk(v,x,scl,smc,method)
% P_REPLICATE_PRK replicates the pre-kernel solution x as a pre-kernel of
% the game space v_sp.
%
% Usage: RepSol=p_replicate_prk(v,x,scl,smc)
% or 
% Usage: [RepSol RepBMat]=p_replicate_prk(v,x,scl,smc)
% This call needs for a 14-person game more than 6 GB disk space.
% One should have at least 10 GB physical memory. Note, adding
% a player increases the memory and disk space requirements by
% a factor of at least four!!!
%
% Define variables:
% output:
%  V_PrkQ                -- Indicates whether the games under consideration
%                           replicate x as a pre-kernel element. 
%  V_SPC                 -- Game space spanned by the basis of the null
%                           space of matrix W.
%  SBCQ                  -- Indicates whether the equivalence class has
%                           changed or not. 
%  SBC                   -- Indicates the set of equivalence class/most
%                           effective coalitions w.r.t. the pre-kernel
%                           element x.
%  Mat_W                 -- Matrix given by equation (7.16) Meinhardt, 2013.
%  P_Basis               -- Basis of the parameter space (null space of mat_W).
%  VarP_Basis            -- Variation in the parameter basis. 
%
%
%  input:
%  v                     -- A Tu-Game v of length 2^n-1. 
%  x                     -- pre-kernel payoff of length(1,n)
%  scl                   -- scaling factor
%  smc                   -- selecting from effc the smallest/largest 
%                           cardinality (optional). Value 1 or 0.
%  method                -- Permissible methods are:
%                           'sparse'  or 
%                            'full'
%                           to compute spares or dense matrices. Default is full.
%                           If you have a memory problem use 'sparse' instead.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   05/18/2011        0.1 alpha        hme
%   10/27/2012        0.3              hme
%   05/12/2014        0.5              hme
%                

if nargin<1
   error('At least the game must be defined!')
elseif nargin==1
   x=p_PreKernel(v);
   scl=1;
   smc=1;
   tol=10^6*eps;
   method = 'full';
elseif nargin==2
   scl=1;
   smc=1;
   tol=10^6*eps;
   method = 'full';
elseif nargin==3 
   smc=1;
   tol=10^6*eps;
   method = 'full';
elseif nargin==4
   if smc > 1, smc=1;
   elseif smc < 0, smc=0;
   else smc=round(smc);
   end
   tol=10^6*eps;
   method = 'full';
else
   if smc > 1, smc=1; 
   elseif smc < 0, smc=0;
   else smc=round(smc);  
   end
   tol=10^6*eps;
end

if isempty(x)
     x=p_PreKernel(v);
else
     prkQ=PrekernelQ(v,x);
   if prkQ==0 
       warning('Input vector is not a pre-kernel element!')
       warning('Input vector must be a pre-kernel element!')
       warning(['Derived game space must be related to the pre-kernel to ' ...
                'have any meaning!'])
   else
   end
end
    

if  nargout>1
  [v_spc, mat_uc, mat_W, ~, A_v, mat_vz]=p_game_space(v,x,scl,smc,method);
  mat_hd=mat_uc';
else
  [v_spc, A_v]=p_game_space_red(v,x,scl,smc,method);
end

sz_v=size(v_spc);
lgsp=sz_v(1);
cA=cell(lgsp,1);
v_spc_prkQ=false(1,lgsp);
sbcQ=false(1,lgsp);

parfor k=1:lgsp
  [~, cA{k}, smat]=BestCoalitions(v_spc(k,:),x,smc);
  lms=abs(smat-smat')<tol;
  v_spc_prkQ(k)=all(all(lms));
  sbcQ(k)=all(all(cA{k}==A_v));
end
cA{lgsp+1}=A_v;


% Formatting output
fields=cell(lgsp+1,1);
parfor k=1:lgsp+1
 fields{k}=strcat('Eq', num2str(k));
end

eqc=cell2struct(cA,fields,1);
if nargout>1
  RepSol=struct('V_PrkQ',v_spc_prkQ,'SBCQ',sbcQ);
  RepBMat=struct('V_SPC', v_spc,'SBC',eqc,'Mat_W',mat_W,'P_Basis', mat_hd,'VarP_Basis', mat_vz);
else
  RepSol=struct('V_PrkQ',v_spc_prkQ,'SBCQ',sbcQ);
end
