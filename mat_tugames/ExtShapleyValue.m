function esh=ExtShapleyValue(v,cs)
% EXTSHAPLEYVALUE computes the extended Shapley-value of a TU-game v w.r.t.
% to a set of restricted coalitions given by cs.
%
% Source: E. Calvoy and Esther Gutiérrez-López (2015); The value in games with restricted cooperation
%         http://www.erices.es/upload/workingpaper/68_0115.pdf
%
% Usage: esh=ExtShapleyValue(v,cs)
%
% Define variables:
%  output:
%  esh       -- The Shapley-value of a TU-game v.
%
%  input:
%  v        -- A TU-Game of length 2^n-1.
%  cs       -- A coalition structure like 
%              cs = [13    51    74    82    97   102   120   127]
%              for {[1 2 5 6] [1 3 4], [1 6 7], [2 3 6 7], [2 4 7], 
%                   [2 5 7], [4 5 6 7], [1 2 3 4 5 6 7]}.
%              The grand coalition must be given.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   12/23/2015        0.8             hme
%



N=length(v);
[~, n]=log2(N);
if (2^n-1)~=N
    error('Game has not the correct size!');
end
if iscell(cs)
   cs=clToMatlab(cs);
end

esh=zeros(1,n);
int=1-n:1:0;
mat=rem(floor(cs(:)*pow2(int)),2)==1;
clS=(mat*ones(n,1))';
fdv=feasible_dividends(v,cs);
shd=fdv./clS;
%rcs=cs;
%rcs(end)=[];
for k=1:n
    Fi=bitget(cs,k)==1;
    if any(Fi)
       esh(k)=sum(shd(Fi));
    end
end

