function v_t=HMS_Reduced_game(v,x,str,tol)
% HMS_REDUCED_GAME computes from (v,x) all Hart/Mas-Colell reduced 
% games on S at x of game v.
%
% Usage: v_t=HMS_Reduced_game(v,x,str)
%
% Define variables:
%
%  output:
%  v_t{1,:} -- All Hart-MasColell reduced games w.r.t. x.
%  v_t{2,:} -- The corresponding Shapley values of all reduced games.
%  v_t{3,:} -- The corresponding sub-coalitions which define a reduced game.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n). Must be efficient.
%  str      -- A string that defines different Methods. 
%              Permissible methods are: 
%              'PRN' that is, the Hart-MasColell reduced game 
%               in accordance with the pre-nucleolus.
%              'PRK' that is, the Hart-MasColell reduced game 
%               in accordance with pre-kernel solution.
%              'SHAP' that is, Hart-MasColell reduced game 
%               in accordance with the Shapley Value, 
%               which is, the original definition. 
%              'MODIC' that is, the Hart-MasColell reduced game 
%               in accordance with the modiclus.
%              'CORE' that is, the Hart-MasColell reduced game 
%               in accordance with the core.
%              Default is 'SHAP'.
%  tol      -- Tolerance value. Its default value is set to 10^4*eps

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/19/2010        0.1 beta        hme
%   06/17/2012        0.2 beta        hme
%   05/27/2013        0.3             hme
%   02/10/2018        0.9             hme
%   04/01/2020        1.9             hme
%   10/17/2021        1.9.1           hme
%                

if nargin<2
  x=ShapleyValue(v);
  n=length(x);
  str='SHAP';
  tol=10^4*eps;  
elseif nargin<3
  n=length(x);
  str='SHAP';
  tol=10^4*eps;  
elseif nargin<4
  n=length(x);
  tol=10^4*eps;  
else
  n=length(x);	
end


N=length(v);  
v_t=cell(3,N-1);

if strcmp(str,'CORE')
%% checking if x is Lorenz maximal within the core.     
   lmQ=LorenzMaxCoreQ(v,x);
else
   lmQ=false;
end


for S=1:N-1
  [v_t{1,S} v_t{2,S} v_t{3,S}]=hms_red_game(v,x,S,n,str,tol,lmQ);
end

%---------------------------------
function [vt subg_sh T]=hms_red_game(v,x,S,n,str,tol,lmQ)


J=1:n;
lmcS=bitget(S,J)==0;
plcS=J(lmcS);
cSpot=2.^(plcS-1);
cS=cSpot*ones(length(plcS),1);
T=SubSets(S,n);
lgt=length(T);
 

vt=zeros(1,lgt);
TorcS=bitor(T,cS);
subT=cell(1,lgt);
subg=cell(1,lgt);
subg_sh=cell(1,lgt);
lg=cell(1,lgt);
plT=cell(1,lgt);
Tz=cell(1,lgt);
sum_py=cell(1,lgt);

for k=1:lgt
 subT{k}=SubSets(TorcS(k),n);
 subg{k}=v(subT{k});
 if strcmp(str,'SHAP')
   subg_sh{k}=ShapleyValue(subg{k});
 elseif strcmp(str,'PRN')
   if length(subg{k})==1
      subg_sh{k}=subg{k};
   else
      try
        subg_sh{k}=msk_prenucl_mod4(subg{k}); % use a thrid party solver instead! 
      catch
        subg_sh{k}=PreNucl(subg{k}); % use default solver instead!
      end
   end
 elseif strcmp(str,'MODIC')
   if length(subg{k})==1
      subg_sh{k}=subg{k};
   else
      try
        subg_sh{k}=msk_modiclus(subg{k}); % use mosek solver.
      catch
        subg_sh{k}=Modiclus(subg{k}); % use default solver instead! 
      end
   end
 elseif strcmp(str,'PRK')
   pl=logical(bitget(TorcS(k),J));
   y=x(pl);    
   subg_sh{k}=PreKernel(subg{k},y);
 elseif strcmp(str,'CORE')
   pl=logical(bitget(TorcS(k),J));
   y=x(pl);
   try
      crQ=CddCoreQ(subg{k});	   
   catch   
      crQ=coreQ(subg{k});
   end   
   if crQ==1
      bcQ=belongToCoreQ(subg{k},y,'rat',tol);
      if bcQ==1      
         subg_sh{k}=y;
      else
         if convex_gameQ(subg{k})
           subg_sh{k}=FindHMConsElm(subg{k},y,TorcS(k),S,n,lmQ); % try to get HM-consistent core element.                
         else            
           %% subg_sh{k}=BaryCenter(subg{k}); % if y is not in the core of the sub-game, we select the bary-center of the core.
           crv=CddCoreVertices(subg{k});
           subg_sh{k}=crv(1,:);  % if y is not in the core of the sub-game, we select the first core element listed.
        end	 
      end	      
   else % taking an non-efficient vector; adding 5 to get also such kind of vector when v(N)=0.
      subg_sh{k}=(ones(1,nnz(pl))*v(N))+5; % not applicable for totally balanced games. 	   
   end	   
 else
   subg_sh{k}=ShapleyValue(subg{k});
 end
 Qk=TorcS(k);
 it=0:-1:1-n;
 lg{k}=rem(floor(Qk(:)*pow2(it)),2)==1;
 plT{k}=J(lg{k});
 Tz{k}=ismember(plT{k},plcS);
 sum_py{k}=Tz{k}*subg_sh{k}';
 vt(k)=v(TorcS(k))-sum_py{k};
end


%%--------Trying to Find Core Element for HM-consistency --------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cr=FindHMConsElm(vQ,xQ,Q,S,n,lmQ)
% FINDHMCONSELM tries to get from (vQ,xQ) a core element of the subgame vQ that is
% monotone w.r.t. xQ.
%
% Source: Dietzenbacher and Sudhoelter, Hart–Mas-Colell consistency and the core in convex games, 
%         IJGT pp. 1-17, (2021). See proof of Lemma 4.1.
%
% Define variables:
%
%  output:
%  cr       -- monotone core element of sub-game vQ w.r.t. xQ,
%              either we have cr <= xQ on all coordinates or at least on
%              the restriciton of R = S \cap Q, i.e., cr_R <= xQ_R.
%
%  input:
%
%  vQ        -- A Tu-Game sub-game vQ of length 2^q-1. 
%  xQ        -- Payoff vector of size(1,q), i.e., restriction of x on Q. 
%               Must be efficient.
%  Q         -- Coalition Q:=T \cup S^c.
%  S         -- Two Person coalition.
%  n         -- Cardinality of the whole player set.
%  lmQ       -- Indicating whether x is Lorenz maximal within the core.
%
%

lq=length(xQ);
try
      if lmQ==1;
         LD_vQ=LorenzSol(vQ);
         crvm=LD_vQ.Cp;
      else
         crvm=CddCoreVertices(vQ);
      end
catch
  crvm=CoreVertices(vQ);
end
greq=crvm<=xQ;
tsm=sum(greq,2);
[lrQ,idxQ]=max(tsm);
if lrQ==lq
   lrg=find(tsm==lrQ);
   scr=crvm(lrg,:);
   pl=1:n;
   pls=pl(logical(bitget(S,pl)));
   plq=pl(logical(bitget(Q,pl)));
   ps=ismember(plq,pls);
   scv=scr(:,ps);
   [~,sidx]=min(sum(scv,2));
   cr=scr(sidx,:);	
%   cr=crvm(idxQ,:);
else %% Let R = S \cap Q, selecting first core element satisfying cr_R <= xQ_R.  
    pl=1:n;
    pls=pl(logical(bitget(S,pl)));
    plq=pl(logical(bitget(Q,pl)));
    ps=ismember(plq,pls);
    rs=xQ(ps);
    ls=length(rs);
    rsl=crvm(:,ps)<=rs;
    srw=sum(rsl,2);
    [lrS,idxS]=max(srw);
    if lrS==ls
       lrg=find(srw==lrS);
       scr=crvm(lrg,:);
       scv=scr(:,ps);
       [~,sidx]=min(sum(scv,2));
       cr=scr(sidx,:);	    
%       cr=crvm(idxS,:);
    else %% exception handling 
       scr=find(srw==lrS); %% getting convex combination.
       lqs=length(scr);
       if lqs > 1
          cr=sum(crvm(scr,:))/lqs;
       else
         cr=crvm(scr,:);
       end
    end
end

