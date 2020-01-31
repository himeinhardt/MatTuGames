function [e,A,smat]=p_BestCoalitions(v,x,smc)
% P_BESTCOALITIONS computes from (v,x,smc) the corresponding
% excess vector, the set of most effective coalitions
% of smallest or largest cardinality depending on the 
% logical value 'smc'. This value must be 0 or 1.
%
%
% Usage: [e A smat]=p_BestCoalitions(v,x,smc)
% Define variables:
% output:
%  e        -- Excess vector
%  A        -- Matrix of most effective coalitions (effc). 
%  smat     -- Matrix of maximal surpluses.
%  bestCell -- Cell contains all sets containing player i but not j. 
%  slcCell  -- Cell contains the set of most effective coalitions. 
%  binCell  -- Cell contains the information of slcCell in bin format.
%
%  input:
%  v      -- A Tu-Game v of length 2^n-1. 
%  x      -- payoff vector of length(1,n)
%  smc    -- selecting from effc the smallest/largest 
%            cardinality (optional). Value 1 or 0.

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   05/18/2011        0.1 alpha       hme
%   09/12/2012        0.2             hme
%   10/22/2012        0.3             hme
%   05/15/2014        0.5             hme
%                


if nargin < 3
   smc=1;
else
   if smc > 1, smc=1; 
   elseif smc < 0, smc=0;
   else  
   end
end

tol=5000*eps;
n=length(x);
% Borrowed from J. Derks
Xm{1}=x(1); for ii=2:n, Xm{1}=[Xm{1} x(ii) Xm{1}+x(ii)]; end
% the excesses of x wrt. the game v
e=v-Xm{1};
clear v Xm;
% Truncate data arrays. 
[se, sC]=sort(e,'descend');
B=eye(n);
smat=-inf(n);
q0=n^2-n;
q=0;
k=1;
pl=1:n;
while q~=q0
  kS=sC(k);
  ai=bitget(kS,pl)==1;
  bj=ai==0;
  pli=pl(ai);
  plj=pl(bj);
  if isempty(plj)==0
    for i=1:numel(pli)
      for j=1:numel(plj)
        if B(pli(i),plj(j))==0
           B(pli(i),plj(j))=k;
           smat(pli(i),plj(j))=se(k); % max surplus of i against j.
           q=q+1;
        end
      end
    end
  end
  k=k+1;
end
m=max(B(:));
e1=se(m)-tol;
le=se>=e1;
tS=sC(le);
te=se(le);
clear se sC;

% Computing the set of most effective coalitions
% and matrix of maximum surpluses.
A=eye(n);
slcCell=cell(n);

% Due to floating point computation the set of most effective
% coalitions might be too small. We set a tolerance value
% hopefully high enough to select a larger set. In case
% that the set of most effective coalitions is not 
% selected correctly, pathological cycles may appear.

parfor i=1:n
   a=bitget(tS,i)==1;
   for j=1:n
     if i<j
       b=bitget(tS,j)==0;
       lij=a & b;
       c_ij=tS(lij);
       ex_ij=te(lij);
       abest_ij=abs(smat(i,j)-ex_ij)<tol;
       slcCell{i,j}=c_ij(abest_ij);
      elseif i>j
       b=bitget(tS,j)==0;
       lij=a & b;
       c_ij=tS(lij);
       ex_ij=te(lij);
       abest_ij=abs(smat(i,j)-ex_ij)<tol;
       slcCell{i,j}=c_ij(abest_ij);
      else
    end
   end
end

smat=tril(smat,-1)+triu(smat,1);

% Assigning the set of selected coalitions to
% matrix A.
parfor i=1:n
  for j=1:n
   if A(i,j)== 0
     lC=length(slcCell{i,j});
     if lC==1
        A(i,j)=slcCell{i,j}; 
     else
         binCell_ij=SortSets(slcCell{i,j},n,lC,smc);
      if smc==1
           A(i,j)=binCell_ij(1);  % Selecting smallest cardinality
       elseif smc==0
           A(i,j)=binCell_ij(end); % Selecting largest cardinality
       else 
           A(i,j)=binCell_ij(1);             % Selecting default
      end
     end
    end
  end
end


%-------------------------------
function Seff=SortSets(effij,n,bd,smc);
% Sorting the set of most effective
% coalitions with respect to their
% cardinality. Ascent ordering.
% Smallest coalitions are coming first.
  Pm=zeros(bd,n);
  for k=1:n, Pm(:,k) = bitget(effij,k);end
  ov=ones(n,1);
  clsize=Pm*ov;
  if smc==1
     mcl=min(clsize);
  else
     mcl=max(clsize);
  end
  eqm=find(clsize==mcl);
  lc=length(eqm);
  if lc~=bd
     effij=effij(eqm);
     Pm=Pm(eqm,:);
     clsize=clsize(eqm);
  end
  pwcl=clsize.^3;
  J=1:n;
  J=J(ones(lc,1),:);
  M=Pm.*J;
  M=M.^(1/2);
  clix=M*ov;
  clnb=clix.*pwcl;
  [~, ix]=sort(clnb);
  Seff=effij(ix);
