function [my_vl, vG]=p_MyersonValue(v,cs,str)
% P_MYERSON_VALUE computes the Myerson value w.r.t. 
% a communication situation using Matlab's PCT.
%
% Usage: my_vl=p_MyersonValue(v,cs,str)
%
% Define variables:
%  output:
%  my_vl    -- The Myerson value of a Tu game.
%  vG       -- The derived restricted game.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
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
%   05/15/2014        0.5             hme
% 
if nargin < 2
   error('A game and communication situation must be given!');
elseif nargin==2
    N=length(v);
    [~, n]=log2(N);
    if (2^n-1)~=N
      error('Game has not the correct size!');
    end
    usQ=p_union_stableQ(cs);
    if usQ==1
       str='us';
    elseif usQ==0
       str='cs';
    else
       str='';
    end
elseif nargin==3
    N=length(v);
    [~, n]=log2(N);
    if (2^n-1)~=N
      error('Game has not the correct size!');
    end
    if strcmp(str,'us')
      usQ=p_union_stableQ(cs);
      if usQ==0
         msg='Coalition structure is not union stable, assuming default.';
         msg2=['Wrong coalition structure will result in an ' ...
              'incorrect result!'];
         warning('usQ:fcs',msg);
         warning('usQ:fcs2',msg2);
         str='';
      end
    elseif strcmp(str,'cs')
      usQ=p_union_stableQ(cs);
      if usQ==1
         msg3='Coalition structure is not a communication structure, assuming union stable.';
         msg4=['Wrong coalition structure will result in an ' ...
              'incorrect result!'];
         warning('usQ:fcs3',msg3);
         warning('usQ:fcs4',msg4);
         str='us';
      end
    else
       usQ=p_union_stableQ(cs);
      if usQ==1
         msg3='Coalition structure is not a communication structure, assuming union stable.';
         msg4=['Wrong coalition structure will result in an ' ...
              'incorrect result!'];
         warning('usQ:fcs3',msg3);
         warning('usQ:fcs4',msg4);
         str='us';
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
