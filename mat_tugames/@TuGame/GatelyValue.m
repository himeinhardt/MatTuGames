function x=GatelyValue(clv,tol)
% GATELYVALUE computes the Gately point of an essential game v. 
%
% Resource: Littlechild and Vaidya (1976) 
%    
% Usage: x=clv.GatelyValue()
% Define variables:
%  output:
%  x        -- Gately Point of the essential game v.
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
%   01/21/2020        1.1             hme
%                

if nargin<2
 tol=10^6*eps;
end
v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;

% upper vector
k=1:n;
vi=clv.tuvi;
Ni=clv.tuSi
uv=v(N)-v(Ni);
%dess=v(N)-sum(vi);
eQ=clv.tuessQ>tol;
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
