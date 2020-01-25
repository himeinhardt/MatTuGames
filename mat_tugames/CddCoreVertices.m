function [core_vert,crst,cr_vol,P]=CddCoreVertices(v,idx,tol)
% CDDCOREVERTICES computes all core vertices of game v, 
% whenever the core exits. The cdd-library by Komei Fukuda is needed.
% It is recommended to install the cdd-library that accompanies
% the Multi-Parametric Toolbox 3.
% http://people.ee.ethz.ch/~mpt/3/
%
% Usage: [core_vert,crst,cr_vol,P]=CddCoreVertices(v)
% Define variables:
%  output:
%  core_vert  -- Matrix of core vertices. Output is numeric.
%  crst       -- The core constraints.
%  cr_vol     -- Core volume, if the core is full dimensional,
%                otherwise zero. 
%  P          -- Returns V- and H-representation (class Polyhedron)
%  input:
%  v          -- A Tu-Game v of length 2^n-1. 
%  idx        -- Specifies the projection plane onto R^3. Default is the empty set.
%                Will be computed internally.
%  tol        -- A positive tolerance value. Its default value is set to 10^9*eps.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   03/07/2011        0.1 beta        hme
%   07/07/2012        0.2 beta        hme
%   10/27/2012        0.3             hme
%   04/24/2014        0.5             hme
%   04/16/2018        1.0             hme
%                

% Here we assume that the user has represented the game correctly.
if nargin<2
  idx='';
  tol=10^9*eps;
elseif nargin < 3
  tol=10^9*eps;
end

N=length(v);
[~, n]=log2(N);

if CddCoreQ(v,tol)==0
  error('Core is empty!');
 else
end

S=1:N;
for k=1:n, PlyMat(:,k) = bitget(S,k);end
%
% Defining core constraints.
%

A1=PlyMat(end,:);
PlyMat(end,:)=[];
A2=PlyMat;
B1=v(N);
v(N)=[];
B2=v';
% Defining the H polyhedron
H=struct('A',[A1;-A2],'B',[B1;-B2],'lin',(1:size(B1,1))');

% Calling cddmex
V=cddmex('extreme',H);
core_vert=V.V;
if nargout == 2
   crst=[H.A,H.B];
elseif nargout > 2
  crst=[H.A,H.B];
  cr_v=core_vert;
  if n>3
     if isempty(idx)
        y=range(cr_v);
        [~, idx]=min(y);
        cr_v(:,idx)=[];
     else
        cr_v(:,idx)=[];
     end
  else
     [X1,X2]=ToSimplex(cr_v);
     cr_v=[X1,X2];
  end 
  P=Polyhedron(cr_v); % V-representation on R^n-1.
  cr_vol=volume(P);
end
