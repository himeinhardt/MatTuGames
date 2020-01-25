function bg=basis_game(n)
% BASES_GAMES computes bases games for n-persons.
% 
% Usage: bg=bases_games(n)
%

% Define variables:
%  output:
%  lb       -- The linear basis of an n-person game. 
%
%  input:
%  n        -- number of players involved.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/18/2013        0.4             hme
%                

N=2^n-1;
int=1-n:1:0;
bg=zeros(N);


for k=1:N
   sS=SubSets(k,n);
   ls=length(sS);
   b=zeros(1,ls);
   mat=rem(floor(sS(:)*pow2(int)),2);
   cS=mat*ones(n,1);
   ec=cS(end);
   for jj=1:ls
       b(jj)=(factorial(ec))/(factorial(cS(jj))*factorial(ec-cS(jj)));
   end    
   bg(k,sS)=b.^(-1);
end


