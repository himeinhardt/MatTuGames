function v=landlord(l,t)
% LANDLORD computes a production game arising from l-landlords and t-tenants.
%
% Source: Rosenmüller (1981) or Driessen (1988).
%
% Usage: v=landlord(l,t)
%
% Define variables:
%  output:
%  v        -- A TU-Game of length 2^n-1.
%
%  output:
%  l        -- Number of landlords given by positive integer.
%  t        -- Number of tenants given by positive integer.   
%
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   12/27/2021        1.9.1           hme
%    
    
n=l+t;
N=2^n-1;
% landlords are located on the first l-th positions, the tenants follow. 
pl=1:n;
ld=1:l;
tnt=l+1:n;

v=zeros(1,N);
for sS=1:N
    mb=pl(logical(bitget(sS,pl)==1));
    lS=nnz(ismember(mb,ld));
    if lS==l
        wts=nnz(ismember(mb,tnt));
        v(sS)=wts^2;
    end    
end    