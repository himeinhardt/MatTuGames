function [acore_vert,acrst,acr_vol,P]=CddAntiCoreVertices(v,idx,tol)
% CDDANTICOREVERTICES computes all anti-core vertices of game v, 
% whenever the anti-core exits. The cdd-library by Komei Fukuda is needed.
% http://www.cs.mcgill.ca/~fukuda/download/cdd and the Matlab interface
% to the cdd solver (cddmex) http://control.ee.ethz.ch/~hybrid/cdd.php.
%
% Usage: [acore_vert acrst]=CddAntiCoreVertices(v)
% Define variables:
%  output:
%  acore_vert  -- Matrix of anti-core vertices. Output is numeric.
%  acrst       -- The anti-core constraints.
%  acr_vol     -- Anti-Core volume, if the anti-core is full dimensional,
%                 otherwise zero. 
%  P           -- Returns V- and H-representation (class Polyhedron)
%
%  input:
%  v           -- A Tu-Game v of length 2^n-1. 
%  idx         -- Specifies the projection plane onto R^3. Default is the empty set.
%                 Will be computed internally.
%  tol         -- A positive tolerance value. Its default value is set to 10^9*eps.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   11/22/2012        0.3             hme
%   07/15/2015        0.7             hme
%                

N=length(v);
[~, n]=log2(N);

if nargin<2
  idx='';
  tol=10^9*eps;
elseif nargin < 3
  tol=10^9*eps;
end


if CddAntiCoreQ(v,tol)==0
  error('Anti-Core is empty!');
 else
end

S=1:N;
for k=1:n, PlyMat(:,k) = bitget(S,k);end
%
% Defining anti-core constraints.
%

A1=-PlyMat(end,:);
A2=PlyMat(1:N-1,:);
B1=-v(N);
B2=v(1:N-1)';
% Defining the H polyhedron
H=struct('A',[A1;A2],'B',[B1;B2],'lin',(1:size(B1,1))');


% Calling cddmex
V=cddmex('extreme',H);
acore_vert=V.V;
%acrst=[H.A,H.B];

if nargout == 2
   acrst=[H.A,H.B];
elseif nargout > 2
  acrst=[H.A,H.B];
  acr_v=acore_vert;
  if n>3
     if isempty(idx)
        y=range(acr_v);
        [~, idx]=min(y);
        acr_v(:,idx)=[];
     else
        acr_v(:,idx)=[];
     end
  else
     [X1,X2]=ToSimplex(acr_v);
     acr_v=[X1,X2];
  end 
  P=Polyhedron(acr_v); % V-representation on R^n-1.
  acr_vol=volume(P);
end

