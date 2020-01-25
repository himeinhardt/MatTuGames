function [Mgc,crst,web_vol,Pweb]=CddWeberSet(v,tol)
% CDDWEBERSET computes the Weber set of game v. Projection method is simplex. 
% Requires the Multi-Parametric Toolbox 3 
% http://people.ee.ethz.ch/~mpt/3/
%
% Usage: [core_vert,crst,cr_vol,P]=CddWeberSetSimplex(v)
% Define variables:
%  output:
%  core_vert  -- Matrix of the vertices of the Weber set. Output is numeric.
%  crst       -- The Weber set constraints.
%  web_vol    -- The volume of the Weber set, if the Weber set is full dimensional,
%                otherwise zero. 
%  Pweb       -- Returns V- and H-representation (class Polyhedron)
%  input:
%  v          -- A Tu-Game v of length 2^n-1. 
%  tol        -- A positive tolerance value. Its default value is set to 10^9*eps.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   06/15/2015        0.7             hme
%


if nargin<2
  tol=10^9*eps;
end

N=length(v);
[~, n]=log2(N);

Mgc=AllMarginalContributions(v);
mgv=ToSimplex3d(Mgc);
Pweb=Polyhedron(mgv);
web_vol=volume(Pweb);
crst=Pweb.H;
