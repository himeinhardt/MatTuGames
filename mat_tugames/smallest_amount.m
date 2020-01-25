function sa=smallest_amount(v)
% SMALLEST_AMOUNT computes the smallest amount vector. 
% This is the smallest amount players can contribute to a coalition.
% 
%  Usage: sa=smallest_amount(v);
%
%
% Define variables:
%  output:
%  r        -- Smallest amount vector
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
%   11/13/2014        0.6             hme
%                
    
if nargin<1
    error('The game must be given!');
else
   N=length(v);
   [~, n]=log2(N);
end

S=1:N;
sa=zeros(1,n);

for i=1:n
 a=bitget(S,i)==1;
 Sa=S(a);
 b=Sa-Sa(1);
 Sa(1)=[];
 bg=b>0;
 s=b(bg);
 vni=v(s);
 sa(i)=min(v(Sa)-vni);
end
