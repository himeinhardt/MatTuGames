function pws=PowerSet(A)
% POWERSET computes the powerset of A without the empty set.
%
% Usage: pws=PowerSet(A)
%
% Define variables:
%  output:
%  pws      -- all subsets of A (power set).
%
%  input:
%  A        -- a set A, like A=[2 3 4]. 
% 
% Example:
% A=[2 3 4];
% pws = PowerSet(A);
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   15/06/2013        0.4             hme

n = length(A);
int=1-n:1:0;
N=2^n-1;
pws=cell(1,N);
S=1:N;
ind=rem(floor(S(:)*pow2(int)),2)==1;
for k=1:N, pws{k} = A(ind(k,:));end

