function [Seff,bmg,cn]=sortsets(S,n)
% SORTSETS sorts an array of coalitions (S,n) given by their
% unique integer representations -- as used in Matlab -- 
% to a canonical order, that is, the coalitions will be sorted according to their cardinality
% and those of equal cardinality are ordered lexicographically.     
% Smallest coalitions coming first, then those of second shortest, etc. 
%
% Example:
%  A collection of sets based on their unique integer representation
%   S=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]
%  will be sorted according to a canonical order given by
%   S=[1 2 4 8 3 5 6 9 10 12 7 11 13 14 15]
%   
%
% Define variables:
%  output:
%  Seff    -- Power set/Coalitions sorted with respect to a canonical order.
%  bmg     -- Sorted sets/coalitions given by their incidence matrix. 
%  cn      -- Canonical numbers of coalitions (sorted in descend order).
%
%  input:
%  S       -- Sets/Coalitions sorted with respect to their unique integer
%              representation.
%  n       -- Maximal number of players. 
% 

%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/29/2010        0.1 alpha       hme
%   04/09/2011        0.1 beta        hme
%   10/14/2012        0.3 beta        hme
%   06/03/2022        1.9.1           hme
%

if nargin<1
  error('The set of coalitions S must be given!');
elseif nargin<2
  error('Cardinality n of the whole player set must be given!');
end



pl=1:n;
bd=length(S);
it=0:-1:1-n;
indM=rem(floor(S(:)*pow2(it)),2);
expl=pl.*ones(bd,n);
imat=indM.*expl;
clm=n-imat;
clm=2.^clm;
dln=sum(clm,2)'; %% canonical numbers of coalitions S.
%% Canonical order R lex T iff d(R) > d(T).
[cn,sid]=sort(dln,'descend');
%% Determining canonical order of coalitions S.
Seff=S(sid);
if nargout == 2
  bmg=imat(sid,:);
elseif nargout == 3
  bmg=imat(sid,:);    	
end  
