function symQ=symmetricQ(clv,tol)
% SYMMETRIC checks if the game v is symmetric.
%
% Usage: symQ=clv.symmetricQ(tol)
%
% Define variables:
%  output:
%  symQ       -- Returns true (1) if the game is symmetric,
%                otherwise false (0).   
%
%  input:
%  clv      -- TuGame class object.
%  tol      -- Tolerance value. Its default value is set to 10^6*eps.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   02/06/2015        0.6             hme
% 

if nargin < 2
   tol=10^6*eps; 
end    

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
it=0:-1:1-n;
S=1:N;
mat=rem(floor(S(:)*pow2(it)),2);
ov=ones(n,1);
csz=mat*ov;

for k=1:n-1
    sS=S(csz==k);
    eQ(k)=all(abs(v(sS)-v(sS(1)))<tol);
end    

symQ=all(eQ);
