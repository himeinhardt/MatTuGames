function [u_coord, sutm]=unanimity_games(v)
% UNANIMITY_GAMES computes the unanimity coordinates or Harsanyi dividends.
% For n>14 this function needs some time to complete.
%
% Usage: [u_coord utmat]=unanimity_games(v)
% Define variables:
%  output:
%  u_coord  -- Unanimity coordinates or Harsanyi dividends.
%  sutm     -- Game basis given as sparse matrix.
%              For a non-sparse matrix, use game_basis() instead.
%
%  input:
%  v        -- A TU-game of length 2^n-1.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/09/2010        0.1 beta        hme
%   10/27/2012        0.3             hme
%   09/30/2013        0.4             hme
%   05/12/2014        0.5             hme
%                



N=length(v);
[~, n]=log2(N);
utmat=false(N);
S=1:N;

for k=1:N
   sS=Subsets(k,n);
   utmat(k,:)=ismembc(S,sS);
end

sutm=sparse(utmat);
if islogical(v)
   v=double(v);
end
u_coord=(sutm\v')';

%--------------------------------------
function sS=Subsets(S,n)

it=0:-1:1-n;
slcP=rem(floor(S(:)*pow2(it)),2)==0;

J=1:n;
sP=J(slcP);

S1=1:S; 

if (2^n-1)==S
  sS=S1;
else
 lsP=length(sP);
 Tni=cell(lsP); 
 for k=1:lsP
  Tni{k}=bitget(S1,sP(k))==0;
 end

 cls=size(Tni);
 ls1=length(S1);
 R=true(1,ls1);
 for k=1:cls(:,2)
  R=Tni{k} & R;
 end
 sS=S1(R);
end
