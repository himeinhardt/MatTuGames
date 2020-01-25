function [vp vp_vec]=veto_players(sv)
% VETO_PLAYERS returns a list of veto players for the simple game sv.
%
% Usage: [vp vp_vec]=veto_players(sv)
% Define variables:
%  output:
%  vp       -- The list of veto players for the simple game sv. It returns
%              a list of zeros whenever the game has no veto player.
%  vp_vec   -- The list of veto players in logical format.
%
%  input:
%  sv       -- A simple game of length 2^n-1.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/15/2010        0.1 beta        hme
%   05/18/2012        0.2 beta        hme
%                

N=length(sv);
[~, n]=log2(N);
vp=zeros(1,n);
Nni=zeros(n,1);

k=1:n;
Nni=bitset(N,k,0);
sp=sv(Nni)==0;
if all(sp==0)
else
 vp(k(sp))=k(sp);
end
vp_vec=vp>0;
