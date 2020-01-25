function [e, A, smat]=BestCoalitions(clv,x,smc)
% BESTCOALITIONS computes from (v,x,smc) the corresponding
% excess vector, the set of most effective coalitions
% of smallest or largest cardinality depending on the 
% logical value 'smc'. This value must be 0 or 1.
%
%
% Usage: [e A smat]=BestCoalitions(clv,x,smc)
% Define variables:
% output:
%  e        -- Excess vector
%  A        -- Matrix of most effective coalitions (effc). 
%  smat     -- Matrix of maximal surpluses.
%
%  input:
%  clv      -- TuGame class object.
%  x        -- payoff vector of length(1,n)
%  smc      -- selecting from effc the smallest/largest 
%              cardinality (optional). Value 1 or 0.

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   10/28/2012        0.3             hme
%                


if nargin < 2
   stx=clv.tustpt;
   if isempty(stx)
     x=clv.select_starting_pt;
   else
     x=stx;
   end
   smc=1;
elseif nargin < 3
   smc=1;
else
   if smc > 1, smc=1; 
   elseif smc < 0, smc=0;
   else  
   end
end

v=clv.tuvalues;
n=clv.tuplayers;
tol=5000*eps;
% the excesses of x wrt. the game v
% Borrowed from J. Derks
Xm{1}=x(1); for ii=2:n, Xm{1}=[Xm{1} x(ii) Xm{1}+x(ii)]; end
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
lcl=length(tS);
te=se(le);
clear se sC;

A=eye(n);
a=false(lcl,n);
c=cell(n);
abest=cell(n);
slcCell=cell(n);
binCell=cell(n);

for i=1:n
   a(:,i)=bitget(tS,i)==1;
end
b=a==0;

% Due to floating point computation the set of most effective
% coalitions might be too small. We set a tolerance value
% hopefully high enough to select a larger set. In case
% that the set of most effective coalitions is not 
% selected correctly, pathological cycles may appear.

for i=1:n-1
   for j=i+1:n
       lij=a(:,i) & b(:,j);
       lji=a(:,j) & b(:,i);
       c{i,j}=tS(lij);
       c{j,i}=tS(lji);
       ex_ij=te(lij);
       ex_ji=te(lji);
       abest{i,j}=abs(smat(i,j)-ex_ij)<tol;
       abest{j,i}=abs(smat(j,i)-ex_ji)<tol;
       slcCell{i,j}=c{i,j}(abest{i,j});
       slcCell{j,i}=c{j,i}(abest{j,i});
   end
end

smat=tril(smat,-1)+triu(smat,1);

for i=1:n-1
  for j=i+1:n
      lCi=length(slcCell{i,j});
      lCj=length(slcCell{j,i});
     if lCi==1
        A(i,j)=slcCell{i,j};
     else
         binCell{i,j}=SortSets(slcCell{i,j},n,lCi,smc);
      if smc==1
           A(i,j)=binCell{i,j}(1);  % Selecting smallest cardinality.
       elseif smc==0
           A(i,j)=binCell{i,j}(end); % Selecting largest cardinality.
      else
           A(i,j)=binCell{i,j}(end);   % Selecting largest cardinality.
      end
     end
     if lCj==1
        A(j,i)=slcCell{j,i};
     else
        binCell{j,i}=SortSets(slcCell{j,i},n,lCj,smc);
       if smc==1
           A(j,i)=binCell{j,i}(1);  % Selecting smallest cardinality.
       elseif smc==0
           A(j,i)=binCell{j,i}(end); % Selecting largest cardinality.
       else
           A(j,i)=binCell{j,i}(end);   % Selecting largest cardinality.
       end
     end
  end
end



%-------------------------------
function Seff=SortSets(effij,n,bd,smc)
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
