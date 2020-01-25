function v=interest_game(dep_vec,iva_vec,i_vec,mths)
% INTEREST_GAME computes from an interest problem the corresponding TU game.
% Source: Lemaire 1991
%
%
% Usage: v=interest_game(dep_vec,iva_vec,i_vec,mths)
%
% Define variables:
%  output:
%  v        -- A TU-interest game.
%  input:
%  dep_vec  -- A deposit vector of positive integers (increasing order). If unsorted
%              sorting will be accomplished by the program.
%  iva_vec  -- An investment vector of length n (amount of money in decreasing order).
%              Determines also the player set.
%  i_vec    -- A vector of the annual interest rate in percent (increasing order). 
%              Must have same length as dep_vec.
%  mths     -- The duration of the investment in months.
%
% Example: From Lemaire p. 78
%
% dep_vec=[1000 3000 5000];
% iva_vec=[1800 900 300];
% i_vec=[.0775 .1025 .12];
% mths=3;
% v=interest_game(dep_vec,iva_vec,i_vec,mths)
%
% v =
%
%   46.1250   17.4375   69.1875    5.8125   53.8125   30.7500   90.0000
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   22/05/2015        0.7             hme
%



if nargin < 4
   mths=12/3;
else 
   mths=12/mths;
end

% Sorting the input vectors.
dep_vec=sort(dep_vec,'ascend');
iva_vec=sort(iva_vec,'descend');
i_vec=sort(i_vec,'ascend');
n2=length(dep_vec);
miv=sum(iva_vec);
% Changing the max deposit value if it is less than the
% total sum of investment.
if dep_vec(n2)<miv;
   dep_vec(n2)=miv+iva_vec(n2)/n2;
end
n=length(iva_vec);
N=2^n-1;
S=1:N;
v=zeros(1,N);
for k=1:n, mat(:,k) = bitget(S,k);end

civa=mat*iva_vec';
edep=[0,dep_vec];
for k=1:n2
    lw=edep(k) <= civa;
    up=civa < edep(k+1);
    sS= up & lw;
    v(sS)=civa(sS)*i_vec(k)/mths;
end
