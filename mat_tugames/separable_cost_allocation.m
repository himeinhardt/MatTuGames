function [sca SC NSC]=separable_cost_allocation(c,str);
% SEPARABLE_COST_ALLOCATION computes from a cost vector c the corresponding
% separable cost allocation. Allowable methods are 'ENSC', 'ACA' and 'SCRB'.
%
% Usage: [sca SC NSC]=separable_cost_allocation(c,str)
% Define variables:
%  output:
%  sca      -- A separable cost allocation of the cost game c depending 
%              on the method. The default method is 'ENSC'.
%  SC       -- The separable (marginal) costs.
%  NSC      -- The non-separable cost in the cost game c.
%
%  input:
%  c        -- A cost game of length 2^n-1. A postive data array is expected.
%  str:     -- A string to invoke one of the following methods:
%  'ENSC'   -- The egalitarian non-separable cost method.
%  'ACA'    -- The alternate cost avoided method.
%  'SCRB'   -- The separable costs remaining benefits method.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/12/2010        0.1 beta        hme
%   05/18/2012        0.2 beta        hme
%   10/27/2012        0.3             hme
%   11/14/2014        0.6             hme
%                



if nargin<2
  str='ENSC'
else
end

N=length(c);
[~, n]=log2(N);

k=1:n;
Si=bitset(N,k,0);

SC=c(N)-c(Si);
NSC=c(N)-SC*ones(n,1);


if strcmp(str,'ENSC')
   sca= SC + (1/n)*NSC;
elseif strcmp(str,'ACA')
   sC=bitset(0,k);
   bet=c(sC)-SC;
   sumb=bet*ones(n,1);
   if sumb==0
     alp=c(sC)/(c(sC)*ones(n,1));
    else
     alp=bet/(sumb);
   end
   sca= SC + alp*NSC;
elseif strcmp(str,'SCRB')
   sC=bitset(0,k);
   v=savings_game(c);
   sca1= scrb_solution(v);
   sca = c(sC) - sca1;
else
   sca= SC + (1/n)*NSC
end
