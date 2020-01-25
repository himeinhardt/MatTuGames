function [caQ,csQ]=compromiseAdmissibleQ(v,tol)
% COMPROMISEADMISSIBLEQ checks if the core cover a TU game v is non-empty.
% In addition, it checks if the game is compromise stable, that is,
% the core cover and the core coincide.
%
% Usage: [caQ,csQ]=compromiseAdmissibleQ(v,tol)
% Define variables:
%  output:
%  caQ      -- Returns 1 (true) whenever the core cover exists, 
%              otherwise 0 (false).
%  csQ      .. Returns 1 (true) whenever the core cover coincide 
%              with the core, otherwise 0 (false).
%
%  input:
%  v        -- A TU-game of length 2^n-1.
%  tol      -- A tolerance value. Default is 10^7*eps
%
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/26/2014        0.5             hme
%                


if nargin < 2
  tol=10^7*eps;
end

N=length(v);
[~, n]=log2(N);
caQ=false;
csQ=false;
% compromise admissible 
[uv,mv]=UtopiaPayoff(v);
ov=ones(n,1);
grQ=all(mv-tol<=uv);
lbQ=mv*ov-tol<=v(N);
ubQ=v(N)-tol<=uv*ov;
caQ = grQ & lbQ & ubQ;

% compromise stable
MV=mv(1); for k=2:n,MV=[MV mv(k) MV+mv(k)]; end
UV=uv(1); for k=2:n,UV=[UV uv(k) UV+uv(k)]; end

sm=sum(uv)-UV;
dv=v(N)-sm; 
ev=max(MV,dv);
csQ=all(v-tol<ev);

