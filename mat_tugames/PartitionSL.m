function pS=PartitionSL(S,cs,n)
% PARTITIONSL determines a partition of S w.r.t. 
% a communication situation.
%
% Usage: pS=PartitionSL(S,cs,n)
%
% Define variables:
%  output:
%  pS       -- Partition of S w.r.t cs.
%
%  input:
%  S        -- A coalition S (Integer).
%  cs       -- A communication situation like [3 5 6]
%              for {[1,2],[1 3],[2 3]}. Other structures induces
%              failures.    
%  n        -- Number of players involved (Integer).    
% 

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/21/2013        0.4             hme
% 
J=1:n;
pl=2.^(J-1);
pS=[];
int=0:-1:1-n;
piS=zeros(1,n);
csm=pl(rem(floor(S(:)*pow2(int)),2)==1);
ca=bitand(S,cs);
rc=ca;
for k=1:n
    ica=unique(ca(ca==pl(k)));
    rc(rc==pl(k))=[];
    if isempty(ica)~=1
       piS(k)=ica;
    end
end
piS=piS(piS>0);
lpi=length(piS);
lrc=length(rc);

if lrc > 1
   rc2=rc;
   st=1;
   e2=0;
   for ii=1:lrc-1
       rc2(1:ii)=[];
       lr2=length(rc2);
       e2=e2+lr2;
       ba(st:e2)=bitand(rc(ii),rc2);
       rc2=rc;
       st=e2+1;
   end
else
   if rc < S
      csm2=pl(rem(floor(rc(:)*pow2(int)),2)==1);
      piS=setdiff(csm,csm2);
      lpi=length(piS);
      ba=0;
   else    
      ba=0;
   end
end  

idx=1;
for m=1:lpi
    imQ=(ba==piS(m));
    if any(imQ)
       isp=[];  
    else
       isp(idx)=piS(m);
       idx=idx+1;
    end    
end

if all(piS==0)
   isp=[];      
end    

if all(ba==0)
   if rc~=0
      rc(rc==0)=[];
      if isempty(rc)~=1
         lp=length(piS);
         if lp >= 1
            for p=1:lp
                nisQ1=any(bitor(piS(p),rc)==rc);
                if nisQ1==1
                   piS(p)=0;  
                end 
             end
             piS(piS==0)=[];
             pS=[piS rc];
         else
              pS=rc;
         end 
      end  
   elseif rc==0
      if lpi==1
         pS=csm;
      else    
         pS=piS;
      end
   else
      rc(rc==0)=[];
      if isempty(rc)~=1    
         lp=length(piS);
         if lp >= 1
            for j=1:lp
                nisQ2=any(bitor(piS(j),rc)==rc);
                if nisQ2==1
                   piS(j)=0;  
                end    
            end
            piS(piS==0)=[];
            pS=[piS rc];
         else
            pS=[piS rc];
         end
      else
         pS=csm;  
      end  
   end   
else
    cf=rc(1); 
    for l=2:lrc
        cf=bitor(cf,rc(l)); 
    end
    fid=find(bitor(isp,cf)==S);
    isp(fid)=[];
    if cf~=S
       pS=[cf,S-cf];
    else 
       pS=sort([isp cf]); 
    end
end    
pS=unique(pS);
