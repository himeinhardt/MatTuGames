function asd_vl=apu_SolidarityValue(v,cs)
% APU_SOLIDARITYVALUE computes the solidarity value w.r.t. 
% the coalition structure cs (a prioir unions).    
%
% Usage: sh_cs=ShapleyValueCS(v,cs)
%
% Define variables:
%  output:
%  asd_vl   -- The solidarity value for coalition structure cs.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  cs       -- A coalition structure like [3 4 8]
%              for {[1,2],[3],[4]}    
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/27/2013        0.4             hme
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
sd_cw=SolidarityValue(cw);
asd_vl=zeros(1,n);
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
         sS=SubSets(cs(jj),n); 
         rcs=cs;
         cbd=cls-1;
         sd_av=zeros(cbd,lcs);
         for l=1:lcm
           plj=pj(l);
           rcs(jj)=plj;
           for i=1:cbd
              av=altered_game(v,rcs,plj,sS(i));
              sd_av(i,:)=SolidarityValue(av);
           end
         end
         % The reduced game. 
         v1=[sd_av(:,jj)',sd_cw(jj)];
         asd_vl(idx)=SolidarityValue(v1);
      else    
         asd_vl(idx)=sd_cw(jj);
      end
  end
else
  asd_vl=SolidarityValue(v);
end    
%-----------------------------------------
function av=altered_game(v,rcs,pj,sS)
% Determines an altered game.
%
% Define variables:
%  output:
%  av       -- An altered game of length 2^lg-1.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  rcs      -- A coalition structure like [3 4 8]
%              for {[1,2],[3],[4]}. Must be a partition.    
%  pj       -- player type pj (an integer).
%  sS       -- A sub-coalition of S.
%

pws=PowerSet(rcs);
lg=numel(pws);
av=zeros(1,lg);
 for k=1:lg
      tQ=any(ismember(pws{k},pj));
      lgc=length(pws{k});
      if lgc==1
         av(k)=v(pws{k});
      elseif tQ==1
         pws2=pws{k};
         pws2(pws2==pj)=[]; 
         cf=pws2(1);
         lcf=length(pws2);
         if lcf==1
           av(k)=v(bitor(cf,sS));  
         else 
           up=bitor(pws2,sS);
           cf=up(1);
           for ii=1:lcf
             c2=up(ii);
             cf=bitor(cf,c2);
           end
           av(k)=v(cf); 
         end
         elseif tQ==0 
         cf=pws{k}(1);  
         for ii=2:lgc
             c2=pws{k}(ii);
             cf=bitor(cf,c2);
         end
         av(k)=v(cf);
      end
 end
 lrc=length(rcs);
 J=1:lrc;
 pws3=PowerSet(J);
 % Correcting order of coalitional values. 
 pm2=zeros(1,lg);
 for j=1:lg
     pm2(j)=sum(2.^(pws3{j}-1));
 end
av=av(pm2);

