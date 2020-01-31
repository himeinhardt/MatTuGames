function MCQ=mex_coalitions(v,x,tol)
% MEX_COALITIONS computes the set of coalitions with maximum excesses
% and checks if this collection is separating. Returns a structure element.
%
% Usage: MCQ=mex_coalitions(v,x,tol)
%
% Define variables:
%  output:
%  mC       -- Field variable that contains the collection of
%              coalitions with maximum excesses.
%  scQ      -- Returns true (1), or false (0). 
%  smQ      -- Matrix of to indicate whether i and j are separating.  
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n).
%  tol      -- Tolerance value. Its default value is set to 10^6*eps.
%

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

if nargin<3
   tol=10^6*eps; 
end     

n=length(x);
N=2^n-1;
e=excess(v,x);
e(N)=[];
[e,sC]=sort(e,'descend');
em=e(1);
e1=em-tol;
le=e>=e1;
mC=sC(le);
[scQ,smQ]=separating_collectionQ(mC,n);
MCQ=struct('mC',mC,'scQ',scQ,'smQ',smQ);
