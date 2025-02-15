function [ma,Mmc]=p_min_aspiration(v)
% P_MIN_ASPIRATION computes the minimum aspiration level of players of game v
% using MATLAB's PCT. 
% 
%  Usage: ma=p_min_aspiration(v);
%
%
% Define variables:
%  output:
%  ma       -- Minimum aspiration level vector.
%  Mmc      -- Maximal marginal contribution vector.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/26/2020        1.9             hme
%                
    
if nargin<1
    error('The game must be given!');
else
   N=length(v);
   [~, n]=log2(N);
end

S=1:N;
sa=zeros(1,n);
Mmc=max(p_AllMarginalContributions(v));
int=1-n:1:0;
k=n:-1:1;
x=fliplr(Mmc);

parfor i=1:n
 a=bitget(S,i)==1;
 Sa=S(a);
 ls=length(Sa);
 mat=rem(floor(Sa(:)*pow2(int)),2)==1;
 mat(:,k(i))=zeros(ls,1);
 sm=mat*x';
 vni=sm';
 ma(i)=max(v(Sa)-vni);
end
