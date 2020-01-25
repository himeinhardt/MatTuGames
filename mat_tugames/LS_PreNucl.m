function x=LS_PreNucl(v)
% LS_PRENUCL computes the least square pre-nucleolus of a game.
%
% Usage: x=LS_PreNucl(v)
% Define variables:
%  output:
%  x        -- Least square pre-nucleolus 
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
%   08/10/2013        0.4             hme
%   09/16/2017        0.9             hme
%     
    
    
N=length(v);
[~, n]=log2(N);
S=1:N;
pl=1:n;
v=gather(v);
sV=arrayfun(@(x) sumSV(v,S,x),pl,'UniformOutput',true);
nV=sum(sV);
sN=2^(n-2);
x=(v(N)/n)+(n*sV-nV)/(n*sN);

%-------------------------
function sV=sumSV(v,S,id)
a=bitget(S,id)==1;
sV=sum(v(a));
