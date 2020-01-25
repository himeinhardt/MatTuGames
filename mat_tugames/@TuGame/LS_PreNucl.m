function x=LS_PreNucl(clv)
% LS_PRENUCL computes the least square pre-nucleolus of a game.
%
% Usage: x=LS_PreNucl(clv)
% Define variables:
%  output:
%  x        -- Least square pre-nucleolus 
%
%  input:
%  clv      -- TuGame class object.
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
    
v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;    


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
