function [vp vp_vec]=veto_players(clv)
% VETO_PLAYERS returns a list of veto players for the simple game sv.
%
% Usage: [vp vp_vec]=veto_players(clv)
%
% Define variables:
%  output:
%  vp       -- The list of veto players for the simple game sv. It returns
%              a list of zeros whenever the game has no veto player.
%  vp_vec   -- The list of veto players in logical format.
%
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
%   10/29/2012        0.3             hme
%                

sv=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
gt=clv.tutype;
if strcmp(gt,'sv')
else
  error('Wrong game type!. Game must be a simple game!')
end

vp=zeros(1,n);
Nni=clv.tuSi;
sp=sv(Nni)==0;

if all(sp==0)
else
 vp(k(sp))=k(sp);
end
vp_vec=vp>0;
