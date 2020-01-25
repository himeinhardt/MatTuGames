function wC=winning_coalitions(mW,n)
% WINNING_COALITIONS computes from a pre-defined set of winning 
% coalitions (e.g. minimal winning coalitions) the whole set of 
% winning coalitions.  
%
% Usage: wC=winning_coalitions(mW,n)
%
% Define variables:
%  output:
%  wC       -- The list of winning coalitions. 
%  input:
%  mW       -- The pre-defined list/vector of winning coalitions.
%  n        -- number of players    
%
%
% Examples:
% Let n=4 and w_coal=[10 12 7 11] be a set of winning coalitions. To
% know whether we have additional winning coalitions, we call
% wC=winning_coalitions(w_coal,4);
% which returns wC=
%            7   10   11   12   13   14   15
% The corresponding simple game can be derived while invoking:
% sv=simple_game(wC,4);
%
% However, for n=5 getting
%    wC2=winning_coalitions(w_coal,5);
% 
%    7   10   11   12   13   14   15   23   26   27   28   29   30   31    
% then typing
%    sv2=simple_game(wC2,5);
% to obtain the corresponding simple game.
% 
    
    

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   05/30/2012        0.2 beta        hme
%   05/15/2014        0.5             hme
%                    

lmW=length(mW);
lmW=sort(lmW);
wCs=cell(lmW,1);
N=2^n-1;
S=1:N;
lsp=1;
stv=0;
for k=1:lmW
    sps=unique(bitor(S,mW(k)));
    wCs{k}=S(ismembc(S,sps));
    lsW=stv+length(wCs{k});
    lsp:lsW;
    wCv(1,lsp:lsW)=wCs{k};
    lsp=lsW+1;
    stv=lsW;
end
wCv=unique(wCv);
wC=S(ismembc(S,wCv));

 
