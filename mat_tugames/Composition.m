function v=Composition(varargin)
% COMPOSITION composite at least 2 TU games up to 15 to a new extended game.
%
% Usage: v=Composition(v1,v2,...,v4) or v=Composition(v1,v2,...,v15) 
%
% Define variables:
%  output:    
%  v         -- A composite TU-game of the inputted games.
%
%  input:
%  vk        -- A collection of TU-games of different size.
%

    
%    
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   30/09/2022        1.9.1           hme
%
     
narginchk(2,15);
n=0;

for k=1:nargin
    Nv(k)=length(varargin{k});
    [~, nv(k)]=log2(Nv(k));
    n=nv(k)+n;
end    
N=2^n-1;
sS=1:N;

if n > 35
   error('Merged markets are too large! Only 35 Agents are allowed in total!');
   return; 
end

for k=1:nargin
    if k==1
       aNv=bitand(sS,Nv(k));
       aNv(aNv==0)=Nv(k)+1;
       gv=[varargin{k},0];
       ints(k,:)=gv(aNv);
    else
       m=sum(nv(1:k-1))+1;
       m2=sum(nv(1:k));
       ply=m:m2;
       M=sum(2.^(ply-1));
       pws=PowerSet(ply);
       sts=clToMatlab(pws);
       uu(sts)=varargin{k};
       aNv=bitand(sS,M);
       aNv(aNv==0)=M+1;
       gv=[uu,0];
       ints(k,:)=gv(aNv);        
    end    
end    
    
    
v=ints(1,:);
for k=2:nargin
    v=v + ints(k,:);
end        
