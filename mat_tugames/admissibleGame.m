function av=admissibleGame(vlm,n)
% ADMISSIBLEGAME computes a symmetric compromise admissible game.
%
% Usage: av=admissibleGame(vlm,n)
% Define variables:
%  output:
%  av       -- A symmetric compromise admissible TU game.
%
%  input:
%  vlm      -- A vector of worths of length n to be assigned to
%              coalitions of the same size.  
%  n        -- Number of players involved (integer)
%
%  Example:
%  Define the vector vlm by
%  vlm=[0 7 12 22] 
%     for 
%     v(S)=0 if |S|=1, v(S)=7 if |S|=2,
%     v(S)=12 if |S|=1, v(S)=22 if |S|=4.
%  Finally, setting n=4, and then invoking 
%     av=admissibleGame(vlm,n)
%  getting
%     av=[0    0    7    0    7    7   12    0    7    7   12    7   12   12   22]
%
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/27/2014        0.5             hme
%                

N=2^n-1;
int=1-n:1:0;
S=1:N;
mat=rem(floor(S(:)*pow2(int)),2)==1; 
clS=mat*ones(n,1);
av=zeros(1,N);
for k=1:n
  cs=S(clS==k);
  av(cs)=vlm(k);
end

