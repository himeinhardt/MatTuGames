function [my_vl, vG]=p_MyersonValue(clv,cs,str)
% P_MYERSON_VALUE computes the Myerson value w.r.t. 
% a communication situation using Matlab's PCT.
%
% Usage: my_vl=p_MyersonValue(clv,cs,str)
%
% Define variables:
%  output:
%  my_vl    -- The Myerson value of a Tu game.
%  vG       -- The derived restricted game.
%
%  input:
%  clv      -- TuGame class object.
%  cs       -- A communication situation like [3 5 6]
%              for {[1,2],[1 3],[2 3]}. Or a union stable
%              structure like [1 7 14 15] for {[1],[1 2 3],
%              [2 3 4],[1 2 3 4]}.
%  str      -- A string that defines different communication structures
%              'us' -- for a union stable structure.
%              ''   -- the empty set otherwise, default.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/07/2013        0.4             hme
%   05/17/2014        0.5             hme
% 

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
if isa(clv,'p_TuVal')
   cms = clv.tu_cs;
   str='cs';
elseif isa(clv,'TuVal')
   cms = clv.tu_cs;
   str='cs';
else
   cms='';
   str='';
end

if nargin < 2
   if isempty(cms)
       error('A game and coalition structure must be given!');
   else
      cs = cms;
      usQ=union_stableQ(cs)
      if usQ==1 
        str='us';
      else
        str='cs';
      end
   end
elseif nargin==2
    usQ=p_union_stableQ(cs);
    if usQ==1
       str='us';
    elseif usQ==0
       str='cs';
    else
       str='';
    end
else
   if strcmp(str,'cs')
      if isempty(cs)
         cs = clv.tu_cs;
      else
         usQ=union_stableQ(cs);
         if usQ==1
            str='us';
         elseif usQ==0
            str='cs';
         else
            str='';
         end
      end
   else
      if isempty(cs)
         us = clv.tu_us;
         if isempty(us)
            msg='Coalition structure is not union stable, assuming default.';
            warning('usQ:fcs',msg);
            cs = clv.tu_cs
            str='';
         else
            cs = us;
         end
      else
         usQ=union_stableQ(cs);
         if usQ==1
            str='us';
         elseif usQ==0
            str='cs';
         else
            str='';
         end
      end
   end
end

if strcmp(str,'us')
int=0:-1:1-n;
vG=zeros(N,1);
 parfor k=1:N
     sb=SubSets(k,n);
     cs=sort(cs);
     iS=sb(ismembc(sb,cs));
     csm=rem(floor(iS(:)*pow2(int)),2)==1;
     scl=csm*ones(n,1);
     mc=max(scl);
     SG=iS(scl==mc);
     vG(k)=sum(v(SG));
 end    
else
 vG=zeros(N,1);
 parfor k=1:N
     SG=PartitionSL(k,cs,n);
     SG(SG==0)=[];
     vG(k)=sum(v(SG));
 end 
end
my_vl=p_ShapleyValue(vG);
