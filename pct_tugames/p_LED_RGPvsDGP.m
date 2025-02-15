function RG=p_LED_RGPvsDGP(v,x,tol)
% LED_RGPvsDGP determines the shift of reduced games of the ECC game or reduced game of the dual cover 
% restricted to N w.r.t. the derived games using Matlab's PCT.  
%
% Source:  H. I. Meinhardt. The Modiclus Reconsidered. Technical report, Karlsruhe Institute of Technology (KIT), Karlsruhe, Germany,
%          2018. URL http://dx.doi.org/10.13140/RG.2.2.32651.75043.
%
%          Meinhardt (2018), "Analysis of Cooperative Games with Matlab and Mathematica".
%
%
% Usage: RG=p_LED_RGPvsDGP(v,x,tol) 
%
% Define variables:
%
%  Output Fields
%  ledvsdgp     -- Returns 1 (true) whenever the shift of all proper reduced games of the ECC game 
%                  w.r.t. the associated proper derived games is identical for all proper coalitions. 
%                  Otherwise 0 (false).
%  rgpvsdgp     -- Returns 1 (true) whenever the shift of all proper reduced games of the dual cover of v
%                  restricted to N  w.r.t. the associated proper derived games is identical for all proper coalitions.
%                  Otherwise 0 (false).
%  t_leddgp     -- Stores the shfits of all coalitions for all compared games within ledvsdgp.
%  t_rgpdgp     -- Stores the shfits of all coalitions for all compared games within rgpvsdgp.
%  mncQ         -- Retunrs 1 (true) whenever the payoff mnc_v is the modiclus, otherwise 0 (false).
%  mnc_v        -- Stores the payoff from which the investigation was imposed.
%  t            -- The sum of the maximal excesses from the primal and dual game.
%
%  Input:
%  v            -- A Tu-Game v of length 2^n-1. 
%  x            -- The modiclus of size(1,n).
%  tol          -- Tolerance value. By default, it is set to 10^6*eps.
%                  (optional)
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   06/08/2018        1.0             hme
%


if nargin<2
   N=length(v);
   tol=10^6*eps;
   mnc_v=msk_modiclus(v);
elseif nargin<3
   N=length(v);
   tol=10^6*eps;
   mncQ=modiclusQ(v,x);
   if mncQ==1
      mnc_v=x;
   else
      warning('Sol:Wrn','Input vector is not the modiclus!');
      mnc_v=x;
   end
else
   N=length(v);
   mncQ=modiclusQ(v,x);
   if mncQ==1 
      mnc_v=x;
   else
      warning('Sol:Wrn','Input vector is not the modiclus!');
      mnc_v=x;
   end
end

MMEX=MMExcess(v,mnc_v);
dc_v=p_DualCover(v);
t=MMEX.mev + MMEX.medv;
[RGP RGPC]=p_Reduced_game_propertyQ(dc_v,[mnc_v,mnc_v],'PRN');
[DGP_mnc_v DGPC_mnc_v]=p_Derived_game_propertyQ(v,mnc_v,'MODIC');
[LEDC_mnc_v LEDCGPQ_mnc_v]=p_Ledcons_propertyQ(v,mnc_v,'MODIC');

dfv=cell(N-1,1);
dfv2=cell(N-1,1);
leddgp=false(1,N-1);
rgpdgp=false(1,N-1);

parfor k=1:N-1
    t1=LEDCGPQ_mnc_v{2}{1,k}-DGPC_mnc_v{2}{k};
    dfv{k}=t1;
    dN1=t1(end);
    t1(end)=[];
    t2=RGPC{2}{1,k} - DGPC_mnc_v{2}{k};
    dfv2{k}=t2;
    dN2=t2(end);
    t2(end)=[];
    if isempty(t1)
       leddgp(k)=abs(dN1)<=tol;
    else
      leddgp(k)=all(abs(t1-2*t)<=tol);
   end
    if isempty(t2)
       rgpdgp(k)=abs(dN2)<=tol;
    else
       rgpdgp(k)=all(abs(t2-t)<=tol);
   end
end

RG.ledvsdgp=leddgp;
RG.rgpvsdgp=rgpdgp;
RG.t_leddgp=dfv;
RG.t_rgpdgp=dfv2;
RG.mncQ=mncQ;
RG.x=mnc_v;
RG.t=t;
