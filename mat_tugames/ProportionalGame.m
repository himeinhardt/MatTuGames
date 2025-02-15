function [v,gv]=ProportionalGame(q,co)
% PROPORTIONALGAME computes from (q,co) a proportional game.
%
% Source: https://math.stackexchange.com/questions/4319847/proportional-coalition-function-coalitional-game-with-fair-splitting
%
% Usage: [v,gv]=ProportionalGame(q,co)
% Define variables:
%  output:
%  v        -- A TU-Game of length 2^n-1.
%  gv       -- Auxiliary TU-Game of length 2^n-1
%
%  output:
%  q        -- Vector of weights of size n.
%  co       -- Cost factor (integer value).
%    
%
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   11/30/2021        1.9.1           hme
%

n=length(q);
N=2^n-1;
pl=1:n;
Ci=2.^(pl-1);
gv=zeros(1,N);
v=zeros(1,N);
cv=ones(1,n)*co;
gv(Ci)=cv.*q;

sS=1:N;
for ii=1:n-1
    a=bitget(sS,ii)==1;
    for jj=ii+1:n
        b=bitget(sS,jj)==1;
        c=bitget(sS,jj)==0;
        Twj=sS(a & c);
        Tj=sS(a & b);
        lT=length(Twj);
        for ll=1:lT 
            gv(Tj(ll))=q(jj)*gv(Twj(ll));
        end    
    end
end   
v=gv-co;
