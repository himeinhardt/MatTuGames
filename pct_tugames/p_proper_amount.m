function pa=p_proper_amount(v)
% PROPER_AMOUNT computes the largest amount players can contribute to a
% proper coalition using MATLAB's PCT.
% 
%  Usage: pa=p_proper_amount(v);
%
%
% Define variables:
%  output:
%  pa       -- Proper amount vector, that is, the largest amount 
%              players can contribute to a proper coalition.
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
%   08/27/2020        1.9             hme
%                
    
if nargin<1
    error('The game must be given!');
else
   N=length(v);
   [~, n]=log2(N);
end

S=1:N-1;
pa=zeros(1,n);

parfor i=1:n
 a=bitget(S,i)==1;
 Sa=S(a);
 b=Sa-Sa(1);
 bg=b>0;
 s=b(bg);
 vni=[0, v(s)];
 pa(i)=max(v(Sa)-vni);
end
