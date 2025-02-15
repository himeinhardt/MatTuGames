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
%   05/17/2022        1.9.1           hme    
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
%% Partition of the Player set.
PPly.sums=smpl; % Sums 
PPly.steps=stpl; % Steps 
PPly.npl=pl(npl); % Null Player
PPly.satv=satv;
