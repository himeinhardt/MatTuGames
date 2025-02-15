function [v,cv]=ComposeMarkets(ptp,seed,tol)
% COMPOSEMARKETS composite from a parameter set defining number of players involved 
% in each separated market situation a market game (merged markets).
%
% Usage: [v,cv]=ComposeMarkets(ptp,seed)
%              or 
%        v=ComposeMarkets;    
%        
%        to get a ten-person composite market game from default parameter.
%    
% Define variables:
%  output:    
%  v          -- A composite market game generated from separated markets.
%
%  input:
%  ptp        -- A collection of integer values, each value specifies the size of a separated market.
%                The sum of the integer values is not allowed to be larger than 35.     
%  seed       -- A seed value, default value is 135.
%  tol        -- A numerical tolerance value, default is set to 10^8*eps.
%
%    
  
    
%
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   30/09/2022        1.9.1           hme
%    
    
if nargin<1
   seed=135;
   rng(seed,'v5uniform'); 
   seed2=ceil(200*rand(1,3));
   ptp=[4 3 3]; %% Constructing 10 player market game.
   lg=length(ptp);
   tol=10^8*eps;   
elseif nargin<2
   lg=length(ptp); 
   seed=135;
   rng(seed,'v5uniform'); 
   seed2=ceil(200*rand(1,lg));
   tol=10^8*eps;
elseif nargin<3
   lg=length(ptp);  
   rng(seed,'v5uniform'); 
   seed2=ceil(200*rand(1,lg)); 
   tol=10^8*eps;    
else
   lg=length(ptp); 
   rng(seed,'v5uniform'); 
   seed2=ceil(200*rand(1,lg)); 
end    

if lg > 15
   error('At most 15 markets can be specified!');
   return;
elseif lg < 2   
   error('At least 2 markets must be specified!');
   return;
end    
tpl=sum(ptp);
if tpl > 35
   error('Merged markets are too large! Only 35 Agents are allowed in total!');
   return; 
end    
    
cv=cell(1,lg);
for k=1:lg
    cv{k}=GetMarketGame(ptp(k),seed2(k),tol);
end    


v=Composite(cv);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v=Composite(mkg)
% COMPOSITE composite at least 2 TU games up to 15 to a new extended game.
%
% Define variables:
%  output:    
%  v          -- A composite TU-game of the inputted games.
%
%  input:
%  mkg        -- A collection of TU-games of different size inputted by a MATLAB cell.
%


%
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   30/09/2022        1.9.1           hme
% 

lm=numel(mkg);
n=0;

for k=1:lm
    Nv(k)=length(mkg{k});
    [~, nv(k)]=log2(Nv(k));
    n=nv(k)+n;
end    
N=2^n-1;
sS=1:N;

for k=1:lm
    if k==1
       aNv=bitand(sS,Nv(k));
       aNv(aNv==0)=Nv(k)+1;
       gv=[mkg{k},0];
       ints(k,:)=gv(aNv);
    else
       m=sum(nv(1:k-1))+1;
       m2=sum(nv(1:k));
       ply=m:m2;
       M=sum(2.^(ply-1));
       pws=PowerSet(ply);
       sts=clToMatlab(pws);
       uu(sts)=mkg{k};
       aNv=bitand(sS,M);
       aNv(aNv==0)=M+1;
       gv=[uu,0];
       ints(k,:)=gv(aNv);        
    end    
end    
    
    
v=ints(1,:);
for k=2:lm
    v=v + ints(k,:);
end 