function [scQ,smQ]=separating_collectionQ(clm,n)
% SEPARATING_COLLECTIONQ verifies if the collection clm is
% separating.
%
% Usage: [scQ,smQ]=separating_collectionQ(clm,n)
%
% Define variables:
% output:
%   scQ       -- Returns true (1), or false (0) 
%   smQ       -- Matrix of to indicate whether i and j are separating.
%
%  input:
%   clm       -- A collection of sets/coalitions represented by their
%                unique integers.
%   n         -- Number of players involved (integer).

%
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   02/13/2015        0.6             hme
%

N=2^n-1;
S=1:N;
for i=1:n
   a(:,i)=bitget(S,i)==1;
end
b=a==0;
smQ=eye(n)==1;
for ii=1:n-1
    for jj=ii+1:n
        tij=S(a(:,ii));
        isQ=clm(ismember(clm,tij));
        if isempty(isQ)==0
           tji=S(b(:,jj));
           smQ(ii,jj)=isempty(clm(ismember(clm,tji)))==0;
        end
        lji=S(a(:,jj));
        isQ2=clm(ismember(clm,lji));
        if isempty(isQ2)==0
           lij=S(b(:,ii));
           smQ(jj,ii)=isempty(clm(ismember(clm,lij)))==0;
        end
    end    
end    
scQ=all(all(smQ));
