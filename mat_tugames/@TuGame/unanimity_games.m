function [u_coord, sutm]=unanimity_games(clv)
% UNANIMITY_GAMES computes the unanimity coordinates or Harsanyi dividends.
% For n>14 this function needs some time to complete.
%
% Usage: [u_coord utmat]=unanimity_games(clv)
%
% Define variables:
%  output:
%  u_coord  -- Unanimity coordinates or Harsanyi dividends.
%  sutm     -- Game basis given as sparse matrix.
%              For a non-sparse matrix, use game_basis() instead.
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
%   10/28/2012        0.3             hme
%   09/30/2012        0.4             hme
%                


v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
sS=cell(1,N);
utmat=zeros(N);
S=1:N;

for k=1:N
   sS{k}=Subsets(k,n);
   utmat(k,:)=ismember(S,sS{k});
end

sutm=sparse(utmat);
u_coord=(sutm\v')';


%--------------------------------------
function sS=Subsets(S,n)

it=0:-1:1-n;
vecS=rem(floor(S(:)*pow2(it)),2);

J=1:n;
slcP=vecS==0;
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
