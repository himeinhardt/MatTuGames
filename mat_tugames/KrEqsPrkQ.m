function kpQ=KrEqsPrkQ(v)
% KREQSPRKQ checks whether the kernel and the pre-kernel coincide 
% for TU-game v using MPT3.
%
%  Source: Chih Chang and Chrong-Hsin Lian (2002)
% 
%  Usage: kpQ=KrEqsPrkQ(v)
%
%
% Define variables:
%  output:
%  prkQ     -- Returns 1 (true) whenever the kernel is equal to
%              the pre-kernel, otherwise 0 (false).
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%    

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   12/20/2014        0.6             hme
%


N=length(v);
[~, n]=log2(N);
k=1:n;
vi=v(bitset(0,k));
essQ=sum(vi)<=v(N);

if essQ==1
   zmQ=zero_monotonicQ(v);
   if zmQ==0
      uv=CddUpperSetVertices(v);
      iv=CddImputationVertices(v);
      tuv=uv';
      Pi=Polyhedron(iv);
      lv=size(uv,1);
      kp=false(1,lv);
      for k=1:lv
         kp(k)=Pi.contains( tuv(:,k) );
      end
      kpQ=all(kp);
   else
      kpQ=zmQ;
   end
else
  kpQ=false;
  warning('KrPQ:01','Game is not essential!');
end
