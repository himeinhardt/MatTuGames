function scrb=scrb_solution(v)
%SCRB_SOLUTION computes separable costs-remaining benefits allocation of game v. 
%
%  USAGE: z=scrb_solution(v)
%  output:
%  scrb      -- A separable cost-remaining benefits allocation of the game v
%
%  input:
%  v        -- A TU game of length 2^n-1.
%
%



%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   11/14/2014        0.6             hme
%


N=length(v);
[~, n]=log2(N);

k=1:n;
Si=bitset(N,k,0);

SC=v(N)-v(Si);
NSC=v(N)-SC*ones(n,1);
scrb=v(N)-v(Si)+NSC/n;
