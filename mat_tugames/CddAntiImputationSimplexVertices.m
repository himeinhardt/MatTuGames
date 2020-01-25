function [imp_vert,crst,ip_vol,P]=CddAntiImputationSimplexVertices(v)
% CDDANTIIMPUTATIONSIMPLEXVERTICES computes all vertices of the anti-imputation set of game v, 
% whenever the anti-imputation set is  essential. Projection method is simplex. 
% The cdd-library by Komei Fukuda is needed.
% http://www.cs.mcgill.ca/~fukuda/download/cdd
%
% It is recommended to install the cdd-library that accompanies
% the Multi-Parametric Toolbox 3.
% http://people.ee.ethz.ch/~mpt/3/
%
% Usage: [imp_vert crst]=CddImputationSimplexVertices(v)
% Define variables:
%  output:
%  imp_vert   -- Matrix of anti-imputation vertices. Output is numeric or a string.
%  crst       -- The anti-imputation constraints.
%  ip_vol     -- Volume of the anti imputation set. 
%  P          -- Returns V- and H-representation (class Polyhedron)
%
%  input:
%  v          -- A Tu-Game v of length 2^n-1. 
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/18/2015        0.7             hme
%                


N=length(v);
[~, n]=log2(N);
if (2^n-1)~=N
   error('Game has not the correct size!');
end

J=1:n;
S=bitset(0,J);
Nk=N-2.^(J-1);
nNk=Nk;
sS=S;
vi=v(S);
s_vi=vi*ones(n,1);
if s_vi<=v(N)
  error('Game is inessential!');
end
Nk(end+1)=N;
it=0:-1:1-n;
PlyMat=rem(floor(Nk(:)*pow2(it)),2);
%
% Defining imputation set.
%
A1=PlyMat(end,:);
A2=PlyMat(1:n,:);
B1=v(N);
B2=v(nNk)';
% Defining the H polyhedron
H=struct('A',[A1;A2],'B',[B1;B2],'lin',(1:size(B1,1))');

% Calling cddmex
V=cddmex('extreme',H);
imp_vert=V.V;
if nargout == 2
   crst=[H.A,H.B];
elseif nargout > 2
  crst=[H.A,H.B];
  ip_v=imp_vert;
  ip_v=ToSimplex3d(imp_vert);
  P=Polyhedron(ip_v); % V-representation on R^n-1.
  ip_vol=volume(P);
end

