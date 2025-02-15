function CH=GetPlayersCharacter(th,w_vec,MmW)
% GETPLAYERSCHARACTER determines from the set of players of the weighted majority game the
% characters Step, Sum and Null Player.
%
% Source:  B Peleg, J. Rosenmueller and P. Sudhoelter (1995); The Kernel and Homogeneous Games with Steps; Chapter 13;
%          Essays in Game Theory: In Honor of M. Maschler; Ed. N. Megiddo
%
%
% Usage: CH=GetPlayersCharacter(th,w_vec);
%
% Define field variables:
%  output:
%  steps    -- Returns all players with character step.
%  sums     -- Returns all players with character sum.
%  npl      -- Returns all null players.
%
%  input:
%  th       -- Threshold to pass a bill (positive number).
%  w_vec    -- Vector of weights (descend ordering).
%  MmW      -- minimal winning coalitions in matrix form (optional).
%              Use function min_homogrep to retrieve correct format.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   05/20/2022        1.9.1           hme
%

n=length(w_vec);
pl=1:n;
N=2^n-1;    
if nargin < 3 % mW is data array
[mW, ~, wmg]=minimal_winning(th,w_vec);
   for k=1:n
       M(:,k)=bitget(mW,k)==1;
   end
   lmW=length(mW);
else % matrix form
   M=MmW;
   clear MmW;
   [lmW,~]=size(M);   
   wmg=weighted_majority(th,w_vec);	
end
NPL=NullPlayers(wmg);
npl=NPL.lnpl;
nnpl=pl(npl==1);
lnpl=length(nnpl);
n1=n-lnpl;
expl=pl.*ones(lmW,n);
imat=M.*expl;
smpl=[];
stpl=[];
fnpl=[];
% Selecting characters 
for k=1:n
    slk=imat(imat(:,k)==k,:);
    lck=size(slk);
    m=zeros(1,lck(1));
    for ii=1:lck(1)
        m(ii)=max(slk(ii,:));        
    end
    lS=min(m);
    mC=sum(w_vec(lS+1:n1));
    isnl=any(k==nnpl);
    if isnl==0
       if w_vec(k) <= mC % Character sum
          smpl=[smpl,k]; 
       else % otherwise character step 
          stpl=[stpl,k]; 
       end
    else
       fnpl=[fnpl,k]; % Null player
    end    
end
CH.steps=stpl;
CH.sums=smpl;
CH.npl=fnpl;
