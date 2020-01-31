function sS=SubSets(S,n)
% SUBSETS computes the power set (subsets) of set S. Must be a number not a
% vector!
%
% Example: 
%     sS=SubSets(9,4) returns all subsets of set/coalition 9, which is
%     1 8 9.
%
% Usage: sS=SubSets(S,n)
% Define variables:
%  output:
%  sS       -- Subsets of the unique integer representation of set S.
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
%   08/08/2010        0.1 beta        hme
%   06/20/2012        0.2 beta        hme
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
slcP=rem(floor(S(:)*pow2(it)),2)==0;

J=1:n;
sP=J(slcP);

S1=1:S; 

if (2^n-1)==S
  sS=S1;
else
 lnsp=length(sP);
 Tni=cell(lnsp); 
  for k=1:lnsp
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

