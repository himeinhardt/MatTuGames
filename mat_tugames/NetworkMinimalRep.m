function NMR=NetworkMinimalRep(E,th,w_vec,tol);
% NETWORKMINIMALREP computes from the set of edges, threshold th and
% the weights w_vec the minimal homogeneous representation of a homogeneous
% network weighted majority game.
%
%
% Usage: mwgs=NetworkMinimalRep(E,th,w_vec,tol);
%
% Define field variables:
%  output:
%  mrQ      -- Returns one (true) whenever a minimal homogeneous representation was found, 
%              otherwise zero (false).
%  hrQ      -- Returns one whenever the representation is homogeneous.
%  mwgs     -- Vector of minimal weights.
%  th       -- Quorum of the minimal homogeneous representation.
%  eQ       -- Returns one whenever the both weighed majority games 
%              i.e., (th0,w_vec0) vs. (th,mwgs), are equal.
%  smpl     -- Players with character sum.
%  stpl     -- Players with character step.
%  npl      -- Players with character null-player.
%  E        -- edge matrix.
%  th0      -- Original quorum.
%  wvec0    -- Original weights.
%
%  input:
%  E        -- An edge matrix of size (lx2) or a cell of numel l.
%              The source must be given by 1, and the sink by the
%              number of the player set. However, the edge matrix
%              can also be of size (lx3) then c can be empty.
%  th       -- Threshold/quorum to pass a bill (positive number).
%  w_vec    -- Vector of weights (descend ordering).
%  hrQ      -- Set one (true) for a homogeneous representation, 
%              otherwise zero (false) or empty set []. Suppresses 
%              the costly check of homogeneity. 
%  tol      -- Numerical tolerance.
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   01/20/2011        1.9             hme
%   06/15/2022        1.9.1           hme
%                    

if nargin<4
   tol=10^8*eps;
end

if nargin <2
   error('The threshold and the weights must be given');
elseif nargin<3
 if iscell(E)
   E=cell2mat(E);
   rw=size(E);
% if the source in E is given by zero.
   me=min(E(:,1));
   if me(1)==0
      s=1;
      E=E+s;
      n=max(E(:,2));
   else   
         n=max(E(:,2)); 
   end
 else
  [m,rw]=size(E);
  if m==2
     E=E';
     me=min(E(:,1));
     if me(1)==0
        E0=E(1:2,:);
        s=1;
        n=max(E(:,2))+s;
        E=E0+s;
     else   
       n=max(E(:,2));          
     end
  else
    me=min(E(:,1));
    if me(1)==0
       E0=E(:,1:2);
       s=1; 
       n=max(E(:,2))+s;
       E=E0+s;
    else   
      n=max(E(:,2));          
    end
  end
 end
 pl=1:n;
 E=E';
 w_vec=ismember(pl,unique(E)'); 
 [wmg,mW,~,n]=NetworkMajorityGame(E,th,w_vec);
else
    [wmg,mW,~,n,E]=NetworkMajorityGame(E,th,w_vec);
    E=E';
end
%% reording the weights in a descend order.
[w_vec,widx]=sort(w_vec,'descend');
lmW=length(mW);
M=false(lmW,n);
%%[mW, ~, wmg]=minimal_winning(th,w_vec);
%% Reordering of the inicence matrix of min-win coalitions.
rpls=zeros(1,n);
for k=1:n
     M(:,k)=bitget(mW,widx(k))==1;
end
qtvec=M*w_vec';
qS=qtvec==th;
hrQ=all(qS);
%%
if hrQ==0
   NMR.mrQ=false;	
   NMR.hrQ=hrQ;  	
   NMR.mwgs=inf;
   NMR.th=inf;
   NMR.eQ=false;  %% games are equal?
   NMR.smpl=[];
   NMR.stpl=[];
   NMR.npl=[];
   NMR.E=E;  %% edge matrix.
   NMR.th0=th;    %% input quorum
   NMR.wvec0(widx)=w_vec;  %% input weights	
   NMR.M=inf;    %% incidence matrix of min-win coalitions.
   return;
end	
NPL=NullPlayers(wmg);
npl=sort(NPL.lnpl);
n=length(w_vec);
pl=1:n;
N=2^n-1;
lmW=length(mW);
%% Selecting non null-players.
if isempty(npl)
    nnpl=pl;
else    
    nnpl=pl(npl==0);
end    
%%
lmW=length(mW);
expl=pl.*ones(lmW,n);
imat=M.*expl;
lnn=length(nnpl);
bd=nnpl(lnn);
Slm=LexMax(mW,n); % lex-max min-win coalition.
Mlm=bitget(Slm,pl)==1;
smpl=[];
stpl=[];
satv=cell(lnn,5);
for ii=1:lnn
    slc=nnpl(ii);
    slk=imat(imat(:,slc)==slc,:);
    lck=size(slk);
    lS=zeros(1,lck(1));
    for jj=1:lck(1)
        lS(jj)=max(slk(jj,:));        
    end
    lk=min(lS);    
    wi=w_vec(slc);
    Dii=lk+1:bd;  %% Domain of player slc.
    mCk=sum(w_vec(Dii));
    if wi<=mCk 
       smpl=[smpl,slc];
       satv{ii,1}=wi;         %% satellite game.
       satv{ii,2}=w_vec(Dii); %% satellite game.
       satv{ii,3}=slc;
       satv{ii,4}=Dii;
       satv{ii,5}=lk;
    else
       stpl=[stpl,slc];
       satv{ii,1}=wi;    %% pseudo satellite game.
       satv{ii,2}=0;     %% pseudo satellite game.
       satv{ii,3}=slc;
       satv{ii,4}=Dii;
       satv{ii,5}=lk;
    end
end
% Defining recursively the minimum weights.
gw_vec=zeros(lnn,1);
gw_vec(lnn)=1;
n1=lnn-1;
stpl1=stpl;
for kk=n1:-1:1
    wi=satv{kk,1};
    slc=satv{kk,3};
    lk=satv{kk,5};
    skk=satv{kk,4};
    if lk+1>slc
       lskk=length(skk);
    else
       lskk=0;	    
       skk=lk;	    
    end	    
    if any(ismember(stpl,kk)) %% character step
       if lskk>0
          ovsk=ones(1,lskk);	       
          gw_vec(kk)=ovsk*gw_vec(skk)+1;	
       else
	  gw_vec(kk)=1;     
       end	       
    elseif any(ismember(smpl,kk)) %% character sum
       mwc=minimal_winning(wi,satv{kk,2});
       if length(mwc)>1
          mwc=LexMax(mwc,n1);  %% Taking the lex-max of the min-win coalitions
       end                     %% of the satellite game of sum kk.
       pl2=1:lskk;
       skk=skk(bitget(mwc(1),pl2)==1);  
       szmw=nnz(skk);                   
       ovsk=ones(1,szmw);
       gw_vec(kk)=ovsk*gw_vec(skk);
    end
end	
g_vec=zeros(n,1);
g_vec(nnpl)=gw_vec;
mwg=zeros(1,n);
mwg(widx)=g_vec';
sz=Mlm*g_vec;
%sz=min(M*g_vec);
wmg_mr=NetworkMajorityGame(E,sz,mwg);
eQ=all(abs(wmg_mr-wmg)<tol);
NMR.mrQ=hrQ && eQ;
NMR.hrQ=hrQ;
NMR.mwgs=mwg;
NMR.th=sz;
NMR.eQ=eQ;
NMR.smpl=widx(smpl);
NMR.stpl=widx(stpl);
NMR.npl=pl(NPL.lnpl);
NMR.E=E;  %% original network 
NMR.th0=th;    %% input quorum
MR.wvec0=zeros(1,n);
MR.wvec0(widx)=w_vec;  %% input weight
NMR.M=M;    %% incidence matrix of min-win coalitions.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lm=LexMax(sS,n)
% LexMax determines the lexicographical maximal coalition from sS, an array of integers
% that represents minimal winning coalitions.
% The input sS is first sorted in accordance with the shortest size,
% and within the block of shortest size the lexicographical maximum is selected.
%
% Usage: lm=LexMax(sS,n);
%
% Define variables:
%  output:
%  lm       -- The lex-max coalition of sS.
%
%  input:
%  sS       -- An array (vector) which contains the information about a set of minimal winning coalitions,
%              for instance, the set of minimal winning coalitions.
%  n        -- Number of player involved, must be an integer.
%

pl=1:n;
it=0:-1:1-n;
indM=rem(floor(sS(:)*pow2(it)),2);
bd=length(sS);
ov=ones(n,1);
clsize=indM*ov;
mcl=min(clsize);
eqm=find(clsize==mcl);
lc=length(eqm);
if lc~=bd
   sS=sS(eqm);
   indM=indM(eqm,:);
end
expl=pl.*ones(lc,n);
imat=indM.*expl;
cl=zeros(1,lc);
for kk=1:lc
    clm=pl(imat(kk,:)>0);
    clm=n-clm;
    clm=2.^clm;
    cl(kk)=sum(clm); %% dual numbers of coalitions sS.
end
[cl,sidx]=max(cl);
%% Determining lexicographic maximal of coalitions sS.
lm=sS(sidx);
