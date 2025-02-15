function [ct,mr_wmg,mr_wgs]=center_solution(th,w_vec,MR,tol)
% CENTER_SOLUTION computes the center solution from the minimal
% representation of an homogeneous weighted majority game.
%
% Source:  B. Peleg, J. Rosenmueller and P. Sudhoelter, The kernel of homogeneous games with steps, in: N.
%          Megiddo, ed., Essays in Game Theory in Honor of Michael Maschler (1994) 163-192.
%          Sudhoelter (1996), Star-shapedness of the kernel for homogeneous games.
%
% Usage: [c,mr_wmg]=center_solution(th,w_vec)
%
% Define field variables:
%  output:
%  ct       -- The center solution of an homogeneous weighted majority
%              game, which is the (pre-)kernel of game mr_wmg.
%  mr_wmg   -- Homogeneous weighted majority game.
%  mr_wgs   -- Minimal representation weights.
%
%  input:
%  th       -- Threshold to pass a bill (positive number).
%  w_vec    -- Vector of weights.
%  MR       -- Structure variable of function min_homogrep (optional)
%  tol      -- Numerical tolerance.
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   09/26/2017        0.9             hme
%   01/20/2021        1.9             hme
%   05/22/2022        1.9.1           hme
%                
if nargin<3
   MR=[];
   tol=10^8*eps;
elseif nargin<4
   tol=10^8*eps;
end
n=length(w_vec);
N=2^n-1;
pl=1:n;
[w_vec,idx]=sort(w_vec,'descend');
if isstruct(MR)==0
   MR=min_homogrep(th,w_vec,[],tol);
end   
if MR.mrQ==0
   ct=inf(1,n);
   mr_wmg=[];
   mr_wgs=[];
   warning('MR:Exit1','Game is not minimal homogeneous. No center solution computable!');   
   return;
end
ct=zeros(1,n);
bd=MR.stpl(end);
wmg=weighted_majority(MR.th,MR.mwgs);
[~,vp]=veto_players(wmg);
wpl=winning_players(MR.th,MR.mwgs);
svp=pl(vp==1);
eqnpl=PlyEqFirstStep(MR.th,MR.mwgs,MR.M);
if MR.hrQ==1 && w_vec(eqnpl.fst)==w_vec(eqnpl.lp) && eqnpl.fst < eqnpl.lp &&  eqnpl.lp < bd && isempty(svp)==1 && isempty(wpl)==1
% %%% standard step game %%%% 
	ms=w_vec(1:eqnpl.lp);
        lm=th-sum(w_vec(eqnpl.lp+1:bd));
        wmg_tr=weighted_majority(lm,ms);
        %% may not work to get the minimal homogeneous representation 
        %% of a non-homogeneous game (lm,ms).
        % MR_tr=min_homogrep(lm,ms,1,tol)
	%ct_tr=MR_tr.mwgs./sum(MR_tr.mwgs); 
        %% rather we rely on Thm 5.4, B. Peleg, J. Rosenmueller and P. Sudhoelter (1994). 
        ct_tr=PreKernel(wmg_tr);
        zv=zeros(1,n-eqnpl.lp);
        ct=[ct_tr,zv]; 
elseif MR.hrQ==1 && isempty(svp)==1 && isempty(wpl)==1
% %%% homogeneous standard game. %%%%%
      c=MR.mwgs./sum(MR.mwgs);
      ct(idx)=c;
elseif MR.hrQ==1 && isempty(svp)==0 && isempty(wpl)==1
   lv=length(svp);   
   ct(svp)=ones(1,lv)/lv;	
elseif MR.hrQ==1 && isempty(svp)==1 && isempty(wpl)==0
    lw=max(wpl);
    n1=n-lw;
    SS=1:N;
    a0=true(1,N);
    a1=true(1,N);
    for jj=1:n
       if jj<lw+1 
          a1=bitget(SS,jj)==0 & a1;
       else
          a0=bitget(SS,jj)==0 & a0;
       end   
    end    
    v0=wmg(SS(a0));
    v1=wmg(SS(a1));
    x1=PreKernel(v1);
    [~,W1]=getMinimalWinning(v1);
    for kk=1:n1 
     M1(:,kk)=bitget(W1,kk)==1; 
    end
    alp=min(M1*x1');
    fac=1/(lw*alp+1);
    t1=alp*ones(1,lw);
    ct=fac*[t1,x1];
else
   ct=inf;  	
end   
mr_wgs=zeros(1,n);
mr_wgs(idx)=MR.mwgs;
if nargout ==2
  mr_wmg=weighted_majority(MR.th,mr_wgs);
elseif nargout == 3
  mr_wmg=weighted_majority(MR.th,mr_wgs);
  mr_wgs.th=MR.th;
  mr_wgs.w=mr_wgs;
end
