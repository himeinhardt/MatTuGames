function RERGP=RE_RGP(v,x,str)
%RE_ RGP checks whether an imputation x is reasonable from both sides for all
% reduced games. 
%
% Usage: RERGP=RE_RGP(v,x,str)
% 
% Define variables:
%
%  output: Fields
%  REQ       -- Returns 1 (true) whenever an imputation x is reasonable from both sides for all reduced games,
%               otherwise 0 (false).  
%  reQ       -- Detailed list of reduced games which satisfy REAS from both sides.
%  RE        -- Returns true (1) if the solution x satisfies REAS,
%               for all reduced games, otherwise false (0).
%  RGP       -- Returns 1 (true) whenever the RGP is satisfied, 
%               otherwise 0 (false).
%  vS        -- All Davis-Maschler or Hart-MasColell reduced games on S at x.
%  impVi     -- Returns a vector of restrictions of x on all S.
%  SMQ       -- Checks if REAS_vS from above is smaller than REAS_vN from above for all reduced games.
%  GRQ       -- Checks if REAS_vS from below is greater than REAS_vN from below for all reduced games.
%  SUB_SMQ   -- Gives the detailed list of SMQ.
%  SLB_GRQ   -- Gives the detailed list of GRQ.
%  N_ub      -- REAS from above for game v.
%  N_lb      -- REAS from below for game v.
%  x         -- Input vector.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1.
%  x        -- payoff vector of size(1,n). Must be efficient.
%  str      -- A string that defines different Methods.
%              Permissible methods are:
%              'PRN' that is, the Davis-Maschler reduced game
%               in accordance with the pre-nucleolus.
%              'PRK' that is, the Davis-Maschler reduced game
%               in accordance with pre-kernel solution.
%              'SHAP' that is, Hart-MasColell reduced game
%               in accordance with the Shapley Value.
%              'MODIC' that is, the Davis-Maschler reduced game
%               equivalence in accordance with the modiclus.
%              'MPRK' that is, the Davis-Maschler reduced game
%               equivalence in accordance with the modified pre-kernel.
%              'PMPRK' that is, the Davis-Maschler reduced game
%               equivalence in accordance with the proper modified pre-kernel.
%              'HMS_PK' that is, Hart-MasColell reduced game
%               in accordance with the pre-kernel solution.
%              'HMS_PN' that is, Hart-MasColell reduced game
%               in accordance with the pre-nucleous.
%              Default is 'PRK'.
%  tol      -- Tolerance value. By default, it is set to 10^6*eps.
%              (optional)
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/06/2018        1.0             hme
%

if nargin<2
  x=PreKernel(v);
  n=length(x);
  tol=10^6*eps;
  str='PRK';
elseif nargin<3
  n=length(x);
  tol=10^6*eps;
  str='PRK';
elseif nargin<4
  n=length(x);
  tol=10^6*eps;
else
  n=length(x);
end

N=length(v);
S=1:N;
it=0:-1:1-n;
PlyMat=rem(floor(S(:)*pow2(it)),2)==1;
rgp_reQ=false(1,N);
%impVec=cell(1,N);
rgpq_sol=cell(1,N);
sol=cell(1,N);


[RGP RGPC]=Reduced_game_propertyQ(v,x,str,tol);
vS=RGPC{2};
impVec=RGPC{4};

REAS(N)=REAS_propertyQ(v,x,tol);
rgp_reQ(N)=REAS(N).reasQ;
ub=REAS(N).ub;
lb=REAS(N).lb;

for k=1:N-1
   N_ub=ub(PlyMat(k,:));
   N_lb=lb(PlyMat(k,:));
   REAS(k)=REAS_propertyQ(vS{1,k},impVec{1,k},tol);
   rgp_reQ(k)=REAS(k).reasQ;
   sub=REAS(k).ub<=N_ub+tol;
   slb=REAS(k).lb>=N_lb-tol;
   sub_sm{k}=sub;
   slb_gr{k}=slb;
   SMQ(k)=all(sub);
   GRQ(k)=all(slb);
end


RERGP.REQ=all(rgp_reQ);
RERGP.reQ=rgp_reQ;
RERGP.RE=REAS;
RERGP.RGP=RGP;
RERGP.vS=vS;
RERGP.impV=impVec;
RERGP.SMQ=SMQ;
RERGP.GRQ=GRQ;
RERGP.SUB_SMQ=sub_sm;
RERGP.SLB_GRQ=slb_gr;
RERGP.N_ub=ub;
RERGP.N_lb=lb;
RERGP.x=x;

