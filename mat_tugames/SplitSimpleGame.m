function SP=SplitSimpleGame(th,w_vec,MR,tol)
% SplitSimpleGame splits a simple game with winning players into sub-games with and without winning players.
%
% Source:  B. Peleg, J. Rosenmueller and P. Sudhoelter, The kernel of homogeneous games with steps, in: N.
%          Megiddo, ed., Essays in Game Theory in Honor of Michael Maschler (1994) 163-192.
%          Sudhoelter (1996), Star-shapedness of the kernel for homogeneous games.
%
% Usage: [sv0,sv1]=SplitSimpleGame(th,w_vec)
%
% Define field variables:
%  output:
%  sv0      -- Sub-game with winning players.
%  sv1      -- Sub-game witout winning players.
%  alpha    -- Reduction Lemma (ii) see Sudhoelter (1996 MSS).
%  fac      -- Reduction Lemma (ii) see Sudhoelter (1996 MSS).
%  t1       -- Incidence vector winning players.
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
%   06/24/2022        1.9.1           hme
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
   warning('MR:Exit1','Game is not minimal homogeneous! No split of the game can be determined!');   
   return;
end
ct=zeros(1,n);
bd=MR.stpl(end);
wmg=weighted_majority(MR.th,MR.mwgs);
[~,vp]=veto_players(wmg);
wpl=winning_players(MR.th,MR.mwgs);
svp=pl(vp==1);
eqnpl=PlyEqFirstStep(MR.th,MR.mwgs,MR.M);

if MR.hrQ==1 && isempty(svp)==1 && isempty(wpl)==0
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
   v0=[];
   v1=[];
   alp=[];
   t1=[];
   fac=[];
end

SP.sv0=v0;
SP.sv1=v1;
SP.alpha=alp;
SP.fac=fac;
SP.t1=t1;
