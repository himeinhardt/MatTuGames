function PPly=PartitionPlySet(th,w_vec)
% PartitionPlySet partitions the set of players of the weighted majority game into
% character Sum, Step, and Null-Player.
%
% Source:  J. Rosenmueller, Homogeneous Games: Recursive Structure and Computation,  
%          Mathematics of Operations Research, Vol. 12, No. 2 (May, 1987), pp. 309-33
%
%
% Usage: PC=PlayersCharacter(th,w_vec);
%
% Define field variables:
%  output:
%  sums     -- Returns all players with character sum.
%  steps    -- Returns all players with character step.
%  nlp      -- Returns all players with character null-player.
%  satv     -- Returns the information of the satellite games.
%
%  input:
%  th       -- Threshold to pass a bill (positive number).
%  w_vec    -- Vector of weights (descend ordering).
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   01/10/2021        1.9             hme
%


[mW, ~, wmg]=minimal_winning(th,w_vec);
NPL=NullPlayers(wmg);
npl=NPL.lnpl;
n=length(w_vec);
pl=1:n;
N=2^n-1;
%% Selecting non null-players.
if isempty(npl)
    nnpl=pl;
else	
    nnpl=pl(npl==0);
end    
for k=1:n 
    M(:,k)=bitget(mW,k)==1;
end
lmW=length(mW);
expl=pl.*ones(lmW,n);
imat=M.*expl;
lnn=length(nnpl);
bd=nnpl(lnn);
smpl=[];
stpl=[];
[~,i0]=min(M(1,:));
i0=i0-1;
satv=cell(lnn,5);
for ii=1:lnn
    slc=nnpl(ii);
    ncl=1:lmW;
    idx=M(:,slc)';
    idx=ncl(idx);
    li=length(idx);
    for jj=1:li
        lS(jj)=max(imat(idx(jj),:));
    end
    lk0=min(lS);
    wi=w_vec(slc);
    %% Looking for fellows (symmetries) of ii.
    %% This allows us to reduce the domain of smaller players.
    we=pl(w_vec==wi);
    t1=we(slc<=we);
    lwe=length(we);
    if lwe > 1
       steqQ=we(end)==n;
       if steqQ
          we(end)=[];
       end	       
    end	    
    if lwe > 1 && slc >= i0
        if isempty(t1)==0	
           t1=t1(1);	
           wsQ=i0<=t1;
        else 
           wsQ=false;
	end	
    elseif lwe > 1 && slc < i0
        wsQ=false;
    elseif lwe == 1 && slc >= i0
        wsQ=true;	    
    elseif lwe == 1 && slc < i0
        wsQ=true;
    else	     
       wsQ=false;	    
    end
    %% Defining domain of smaller players.    
    if wsQ==0
       lk=max(lk0,we(end));
    else
       if slc == we(1)  && i0 < we(1) && lwe == 1 && lk0-t1>=2
          lk=lk0;
       elseif slc == we(1)  && i0 < we(1) && lwe == 1 && lk0-t1<2
          lk=we(1);
       elseif slc >= we(1) && i0 <= we(end)  && lwe > 1
           lk=max(lk0,we(end));
       else
          lk=lk0;
       end
    end
    Dii=lk+1:bd;  %% Domain of smaller players.
    mCk=sum(w_vec(Dii));
    if wi<=mCk && slc < lk+1 
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
%% Partition of the Player set.
PPly.sums=smpl; % Sums 
PPly.steps=stpl; % Steps 
PPly.npl=pl(npl); % Null Player
PPly.satv=satv;
