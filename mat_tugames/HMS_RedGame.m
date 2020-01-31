function [vt subg subg_sh]=HMS_RedGame(v,x,S)
% HMS_REDGAME computes from (v,x,S) a Hart-MasColell reduced game vS on S 
% at x for game v.
%
% Usage: [vt subg subg_sh]=HMS_RedGame(v,x,S)
% Define variables:
%  output:
%  vt       -- The Hart-MasColell reduced game vS w.r.t. x.
%  subg     -- The corresponding sub-game of v on S to define 
%              the reduced game vS.
%  subg_sh  -- The list of Shapley values w.r.t. all sub-games subg. 
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n).
%  S        -- A coalition/set identified by its unique integer representation.
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/21/2010        0.1 beta        hme
%   06/19/2012        0.2 beta        hme
%   10/27/2012        0.3             hme
%                


N=length(v);
n=length(x);
J=1:n;
lmcS=bitget(S,J)==0;
plcS=J(lmcS);
cSpot=2.^(plcS-1);
cS=cSpot*ones(length(plcS),1);
T=SubSets(S,n);
lgt=length(T);
 

if S==N
  vt=v;
  subg=0;
  subg_sh=0;
else

 vt=zeros(1,lgt);
 TorcS=bitor(T,cS);

 subT=cell(1,lgt);
 subg=cell(1,lgt);
 subg_sh=cell(1,lgt);
 lg=cell(1,lgt);
 plT=cell(1,lgt);
 Tz=cell(1,lgt);
 sum_py=cell(1,lgt);

 for k=1:lgt
  subT{k}=SubSets(TorcS(k),n);
  subg{k}=v(subT{k});
  subg_sh{k}=ShapleyValue(subg{k});
  it=0:-1:1-n;
  Qk=TorcS(k); 
  lg{k}=rem(floor(Qk(:)*pow2(it)),2)==1;
  plT{k}=J(lg{k});
  Tz{k}=ismember(plT{k},plcS);
  sum_py{k}=Tz{k}*subg_sh{k}';
  vt(k)=v(TorcS(k))-sum_py{k};
 end
end
