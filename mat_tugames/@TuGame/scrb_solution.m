function scrb=scrb_solution(clv)
%SCRB_SOLUTION computes separable costs-remaining benefits allocation of game v. 
%
%  USAGE: z=scrb_solution(v)
%  output:
%  scrb      -- A separable cost-remaining benefits allocation of the game v
%
%  input:
%  clv      -- TuGame class object.
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


v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;

k=1:n;
Si=bitset(N,k,0);

SC=v(N)-v(Si);
NSC=v(N)-SC*ones(n,1);
scrb=SC+NSC/n;
