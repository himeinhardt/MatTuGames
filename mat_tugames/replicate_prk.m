function RepSol=replicate_prk(v,x,scl,smc)
% REPLICATE_PRK replicates the pre-kernel solution x as a pre-kernel of
% the game space v_sp.
%
% Usage: RepSol=replicate_prk(v,x,scl,smc)
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
%  Mat_W                 -- Matrix given by equation (5.15) Meinhardt, 2010.
%  P_Basis               -- Basis of the parameter space (null space of mat_W).
%  VarP_Basis            -- Variation in the parameter basis. 
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
%   02/05/2011        0.1 beta        hme
%   03/02/2011        0.1 beta 2      hme
%                

if nargin<1
   error('At least the game must be defined!')
elseif nargin==1
   x=PreKernel(v);
   scl=1;
   smc=1;
   tol=10^6*eps;
elseif nargin==2
   scl=1;
   smc=1;
   tol=10^6*eps;
elseif nargin==3 
   smc=1;
   tol=10^6*eps;
else
   if smc > 1, smc=1; 
   elseif smc < 0, smc=0;
   else smc=round(smc);  
   end
   tol=10^6*eps;
end

if isempty(x)
     x=PreKernel(v);
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
    
[v_spc mat_uc mat_W mat_V A_v mat_vz]=game_space(v,x,scl,smc);
mat_hd=mat_uc';
sz_v=size(v_spc);
lgsp=sz_v(1);
cA=cell(lgsp,1);
v_spc_prkQ=zeros(1,lgsp);
sbcQ=zeros(1,lgsp);

for k=1:lgsp
   [e cA{k,1} smat]=BestCoalitions(v_spc(k,:),x,smc);
   lms=abs(smat-smat')<tol;
   v_spc_prkQ(k)=all(all(lms));
   sbcQ(k)=all(all(cA{k,1}==A_v));
end
cA{lgsp+1,1}=A_v;
% Formatting output
fields=cell(lgsp+1,1);
for k=1:lgsp+1
 fields{k,1}=strcat('Eq', num2str(k));
end
eqc=cell2struct(cA,fields,1);
RepSol=struct('V_PrkQ',v_spc_prkQ,'V_SPC', v_spc,'SBCQ',sbcQ,'SBC',eqc,'Mat_W',mat_W,'P_Basis', mat_hd,'VarP_Basis', mat_vz);

