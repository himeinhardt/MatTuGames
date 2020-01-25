function r=reasonable_outcome(v)
%REASONABLE_OUTCOME computes the largest amount vector. 
%This is the largest amount players can contribute to a coalition.
% 
%  Usage: r=reasonable_outcome(v);
%
%
% Define variables:
%  output:
%  r        -- Largest amount vector
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
%   05/14/2012        0.2 beta        hme
%   10/12/2012        0.3             hme
%                
    
if nargin<1
    error('The game must be given!');
else
   N=length(v);
   [~, n]=log2(N);
end

S=1:N;
r=zeros(1,n);

for i=1:n
 a=bitget(S,i)==1;
 Sa=S(a);
 b=Sa-Sa(1);
 bg=b>0;
 s=b(bg);
 vni=[0, v(s)];
 r(i)=max(v(Sa)-vni);
end
