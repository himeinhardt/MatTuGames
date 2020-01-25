function pS=PartitionSA(S,hs,n)
% PARTITIONSA determines a partition of S w.r.t. 
% a communication situation.
%
% Usage: pS=PartitionSA(S,hs,n)
%
% Define variables:
%  output:
%  pS       -- Partition of S w.r.t hypergraph communication situation.
%
%  input:
%  S        -- A coalition S (Integer).
%  hs       -- A communication situation like [3 5 6]
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
%   07/31/2013        0.4             hme
% 
J=1:n;
pl=2.^(J-1);
% Separating coalition S w.r.t. hs.
ca=bitand(S,hs);
lc=length(ca);
rc=ca;
if lc == 1
   sp=ca;
   is=[];
elseif lc == 2
   ba=bitand(ca(1),ca(2));
   if ba == 0
      sp=ca;
      is=[];
   else
      is=bitor(ca(1),ca(2));
      sp=[];
   end    
else 
  ugr=lc-1;
  ba=cell(1,ugr);
  sp=cell(1,ugr);
  is=cell(1,ugr);   
  for ii=1:ugr
      rc(1)=[];
      ba{ii}=bitand(ca(ii),rc);
      sp{ii}=rc(ba{ii}==0);
      is{ii}=[ca(ii),rc(ba{ii}>0)];
  end
end

if iscell(is)
  sp=unique(cell2mat(sp));
  is=unique(cell2mat(is));
  lis=length(is);
  glis=lis-1;
  slc=zeros(1,glis);
     for ii=1:glis
         uc=is(ii);
         for jj=ii+1:lis
             uc=bitor(uc,is(jj));
         end
         slc(ii)=uc;
     end
else
   slc=is;    
end  
slc=unique(slc);
ls=length(slc);
if ls > 1
   ugr2=ls-1;
   slc2=zeros(1,ugr2);
   for ii=1:ugr2
      uc2=slc(ii);
      for jj=ii+1:ls 
          ba2=bitand(uc2,slc(jj));
          if ba2 > 0
             uc2=bitor(slc(ii),slc(jj)); 
          end    
      end
      slc2(ii)=uc2;
   end
   slc=slc2;
else   
end
pS=sort([slc sp]);
% Selecting unitary coalitions 
lpS=length(pS);
ui=zeros(1,n);
for k=1:lpS
    idx=J(bitget(pS(k),J)==1);
    ui(idx)=idx;
end
piS=pl(ui==0);
pS=sort(unique([piS,pS]));
