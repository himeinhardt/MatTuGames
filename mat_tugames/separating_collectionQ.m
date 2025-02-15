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
%                unique integers. A row (vector) of positive integers.
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
%   01/23/2022        1.9.1           hme
%

smQ=eye(n)==1;
for ii=1:n-1
    for jj=ii+1:n
        aii=bitget(clm,ii)==1;
        bjj=bitget(clm,jj)==0;
        inj=any(aii & bjj);
        ajj=bitget(clm,jj)==1;
        bii=bitget(clm,ii)==0;
        eqa=all(aii==ajj);
        eqb=all(bii==bjj);
        jni=any(ajj & bii);
        if inj && jni || (eqa && eqb)==1;
           smQ(ii,jj)=true;
           smQ(jj,ii)=true;
        else
           smQ(ii,jj)=false;
           smQ(jj,ii)=false;
        end    
    end    
end    
scQ=all(all(smQ));
