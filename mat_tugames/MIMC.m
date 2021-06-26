function [betv,lcv,A]=MIMC(v,tol)
% MIMC computes the vector of minimum increase in players marginal contribution when they leave the grand coalition. 
%
% Source: Multi-Player Allocations in the Presence of Diminishing Marginal Contributions:
%         Cooperative Game Analysis and Applications in Management Science, 2020, to appear in Management Science 
%         Thm. 5
%
%  Usage: [betv,lcv,A]=MIMC(v,tol)
%
% Define variables:
%  output:
%  betv     -- Vector of minimum increase in players marginal contribution when they leave the grand coalition. 
%  lcv      -- Least core value whenever the game satisfies the ADMC property, otherwise -inf.
%  A        -- A pair matrix, which indicates the increase in player i's marginal contribution when he leaves the grand coalition for a 
%              (n-1) coalition without player j.  
%
%  input:
%  v        -- A TU-game of length 2^n-1.
%  tol      -- Tolerance value. By default, it is set to 10^6*eps.
%              (optional) 


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   05/26/2020        1.9             hme
%                

if nargin<2
   tol=10^6*eps;
end

N=length(v);
[~, n]=log2(N);
k=1:n;
Ni=N-2.^(k-1);
betv=(inf(1,n));
A=diag(betv);

for i=1:n
  for j=1:n
   if A(i,j)==0	  
      Nni=bitset(N,i,0);
      Nnj=bitset(N,j,0);
      Nnij=bitset(Nni,j,0);
      A(i,j)=v(Nnj)-v(Nnij) - v(N)+v(Nni);
   end
 end
end

for i=1:n
  betv(i)=min(A(i,:));
end

if nargout >= 2
 if ADMC_gameQ(v)
% upper payoff/Separable Contribution to the Grand Coalition
  bv=v(N)-v(Ni);
% Non-Separable Contribution/Cost (non-negative/non-positive)
  nsc=v(N)-sum(bv);
% Inequality 9
  grQ=betv>=nsc/(n-1);
  if any(grQ) & ADMC_CoreQ(v,tol)==0
     lcv=nsc/n;
  else
     lcv=0;  	
  end
 else
   lcv=-inf;
 end
end
