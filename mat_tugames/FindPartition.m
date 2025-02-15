function PartQ=FindPartition(v,x,smc,tol)
% FINDPARTITION tries to find a partition and anti partition w.r.t. payoff vector x. 
%
% Usage: PartQ=FindPartition(v,x,smc,tol)
%
% Define variables:
%  output:
%  ptn         -- This structure element returns the set of most effective coalitions 
%                 with largest excess (best coalitions).
%  ptnQ        -- This structure element returns 1 (true) if ptn forms a partition of the player set, 
%                 otherwise 0 (false). 
%  aptn        -- The anti partition if pnt is a partition, otherwise none.
%  setD        -- The set of largest excess other than the grand coalition w.r.t. x.  
%  lex         -- Largest excess other than the grand coalition.
%
%    
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n) (optional), default is the pre-nucleolus.
%  smc      -- The cardinality of the set of most effective coalitions.
%              Smallest and largest cardinality can be chosen.
%              Permissible values are:
%              1 to invoke the smallest cardinality. 
%              0 to invoke the largest cardinality.
%  tol      -- A tolerance value, default is 10^8*eps;
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   12/04/2013        0.5             hme
%   01/25/2022        1.9.1           hme
%


if nargin < 2
  try 
     x=prenucl(v);
  catch
     x=glpk_prenucl(v);
  end
  smc=1;
  tol=10^8*eps;
elseif nargin < 3
  smc=1;
  tol=10^8*eps;
elseif nargin < 4
  tol=10^8*eps;
end


[exc, A, smat]=BestCoalitions(v,x,smc);
N=length(v);
exc(end)=[];
lex=max(exc);
tm=abs(smat-lex)<tol;
sm1=tril(tm,-1);
sm2=triu(tm,1);
fp1=A(sm1)';
fp2=A(sm2)';
lp1=length(fp1);
lp2=length(fp2);
if lp1==lp2
   mf=[fp1;fp2];
else
   mf=[unique(fp1);unique(fp2)];
end    
sNQ=sum(mf,1)==N;
aNQ=any(sNQ);
C=unique(mf(:,sNQ),'rows')';
R=unique(C,'rows');
fp=unique([fp1,fp2]);
NQ=sum(fp)==N;
cls=fp;
cls(1)=[];
ffp=fp(1);
mbaQ=[];
while 1
  if isempty(cls)
     break;
  end
  baQ=bitand(ffp,cls);
  mbaQ=[mbaQ, baQ];
  ffp=cls(1);
  cls(1)=[];
end
aQ=all(mbaQ==false);
fpQ=NQ & aQ;
if aNQ==1
   fp=R';  
   aptn=N-fp;
   fpQ=aNQ;
elseif fpQ==1
   fp=fp'; 	
   aptn=N-fp;
else
   aptn=fp; % Checking if from a supposed anti-partition one gets a partition.
   fp=N-fp;
   NQ=sum(fp)==N;
   cls=fp;
   cls(1)=[];
   ffp=fp(1);
   mbaQ=[];
   while 1
   if isempty(cls)
      break;
   end
   baQ=bitand(ffp,cls);
   mbaQ=[mbaQ, baQ];
   ffp=cls(1);
   cls(1)=[];
   end
   aQ=all(mbaQ==false);
   fpQ=NQ & aQ;
end
setD=SetofLargestExcess(v,x,tol);
PartQ=struct('ptn',fp,'ptnQ',fpQ,'aptn',aptn,'setD',setD,'lex',lex);
