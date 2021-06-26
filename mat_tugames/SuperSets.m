function sS=SuperSets(S,n)
% SUPERSETS computes the super-sets of set S. Must be a number not a
% vector!
%
% Example: 
%     sS=SuperSets(9,4) returns all super-sets of set/coalition 9, which is
%     9   11   13   15   25   27   29   31.
%
% Usage: sS=SuperSets(S,n)
% Define variables:
%  output:
%  sS       -- Super-sets of the unique integer representation of set S.
%
%  input:
%  S        -- A positive number, which represents a unique 
%              integer representation of a set S.
%  n        -- A positive number, which indicates the number of 
%              players in a Tu-game

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   03/28/2021        1.9             hme
%                

if nargin<1
  error('A coalition S and the maximal number of players are missing! Both inputs must be an integer!');
 elseif nargin==1
  if length(S)>1
    error('Coalition must be represented by an integer!');
  else 
    error('The maximal number of players must be given! Must be an integer.');
  end
 else
  if length(S)>1
    error('Coalition must be represented by an integer!');
  elseif length(n)>1 
    error('The maximal number of players must be given by an integer!');
  else
  end
end
 


it=0:-1:1-n;
slcP=rem(floor(S(:)*pow2(it)),2)==1;

J=1:n;
sP=J(slcP);
N=2^n-1;
S1=S:N; 

if N==S
  sS=S1;
else
 lnsp=length(sP);
 Tni=cell(lnsp); 
  for k=1:lnsp
    Tni{k}=bitget(S1,sP(k))==1;
  end

cls=size(Tni);
ls1=length(S1);
R=true(1,ls1);
 for k=1:cls(:,2)
  R=Tni{k} & R;
 end
sS=S1(R);
end

