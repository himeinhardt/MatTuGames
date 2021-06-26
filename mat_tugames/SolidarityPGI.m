function [apgi vr]=SolidarityPGI(v,cs)
% SolidarityPGI computes the solidarity Holler index w.r.t. 
% a priori unions cs.
%
% Usage: apgi=SolidarityPGI(v,cs)
%
% Define variables:
%  output:
%  apgi     -- The a priori unions Holler index.
%
%  input:
%  v        -- A simple Tu-Game v of length 2^n-1. 
%  cs       -- A coalition structure like [3 4 8]
%              for {[1,2],[3],[4]}. Must be a partition.    
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   12/12/2020        1.9             hme
% 
if nargin < 2
   error('A game and coalition structure must be given!'); 
elseif nargin==2
    N=length(v);
    [~, n]=log2(N);
    if (2^n-1)~=N
      error('Game has not the correct size!');
    end
end    

lcs=length(cs);
pws=PowerSet(cs);
% Constructing intermediate game with a priori unions (quotient
% game). 
lg=numel(pws);
if lg > 1
cw=zeros(1,lg);
 for k=1:lg
     lgc=length(pws{k});
      if lgc==1
         cw(k)=v(pws{k});
      else
         cf=pws{k}(1);  
         for ii=2:lgc
             c2=pws{k}(ii);
             cf=bitor(cf,c2);
         end
         cw(k)=v(cf);
      end
 end
pu=1:lcs;
pws2=PowerSet(pu);
pcl2=zeros(1,lg);
% Correcting order of coalitional values. 
for k=1:lg
    pcl2(k)=sum(2.^(pws2{k}-1));
end
cw=cw(pcl2);
sh_cw=PGI_SV(cw);
apgi=zeros(1,n);
J=1:n;
int=0:-1:1-n;
csm=rem(floor(cs(:)*pow2(int)),2)==1;
  for jj=1:lcs
      Tj=csm(jj,:);
      idx=J(Tj);
      pj=2.^(idx-1);
      lcm=length(pj);
      cls=2^lcm-1;
      if lcm > 1
         apgi([idx])=sh_cw(jj)/lcm;
      else    
         apgi([idx])=sh_cw(jj);
      end
  end
else
   apgi=PGI_SV(v); 
end  
