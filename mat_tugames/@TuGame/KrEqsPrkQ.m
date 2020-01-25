function kpQ=KrEqsPrkQ(clv)
% KREQSPRKQ checks whether the kernel and the pre-kernel coincide 
% for TU-game v using MPT3.
%
%  Source: Chih Chand and Chrong-Hsin Lian (2002)
% 
%  Usage: kpQ=clv.KrEqsPrkQ()
%
%
% Define variables:
%  output:
%  prkQ     -- Returns 1 (true) whenever the kernel is equal to
%              the pre-kernel, otherwise 0 (false).
%  input:
%  clv      -- TuGame class object.
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

essQ=clv.tuessQ;

if essQ==1
   zmQ=clv.zero_monotonicQ;
   if zmQ==0
      uv=clv.CddUpperSetVertices;
      iv=clv.CddImputationVertices;
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
