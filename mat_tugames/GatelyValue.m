function x=GatelyValue(v,tol)
% GATELYVALUE computes the Gately point of an essential game v. 
%
% Resource: Littlechild and Vaidya (1976) 
%    
% Usage: x=GatelyValue(v)
% Define variables:
%  output:
%  x        -- Gately Point of the essential game v.
%
%  input:
%  v        -- An essential TU-game of length 2^n-1.
%
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   01/21/2020        1.1             hme
%                

if nargin<2
 tol=10^6*eps;
end
%tol=-tol;
N=length(v);
[~, n]=log2(N);
% upper vector
k=1:n;
vi=v(2.^(k-1));
Ni=bitset(N,k,0);
uv=v(N)-v(Ni);
dess=v(N)-sum(vi);
eQ=dess>tol;
x=-inf(1,n);

 if eQ==1
    d=(sum(uv)-v(N))/dess;
    if d~=-1
       x=(uv+d*vi)./(d+1);
    end
 else
   warning('Gately:Exit','Game is not essential, Gately Point does not exist!');    
 end
end