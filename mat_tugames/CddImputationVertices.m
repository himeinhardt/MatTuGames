function [imp_vert,crst,ip_vol,P]=CddImputationVertices(v,idx)
% CDDIMPUTATIONVERTICES computes all vertices of the imputation set of game v, 
% whenever the imputation set is  essential. The cdd-library by 
% Komei Fukuda is needed.
% http://www.cs.mcgill.ca/~fukuda/download/cdd
%
% It is recommended to install the cdd-library that accompanies
% the Multi-Parametric Toolbox 3.
% http://people.ee.ethz.ch/~mpt/3/
%
% Usage: [imp_vert crst]=CddImputationVertices(v)
% Define variables:
%  output:
%  imp_vert   -- Matrix of imputation vertices. Output is numeric or a string.
%  crst       -- The core constraints.
%  ip_vol     -- Volume of the imputation set. 
%  P          -- Returns V- and H-representation (class Polyhedron)
%
%  input:
%  v          -- A Tu-Game v of length 2^n-1. 
%  idx        -- Specifies the projection plane onto R^3. Default is the empty set.
%                Will be computed internally.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   09/11/2010        0.1 beta        hme
%   09/18/2012        0.2 beta        hme
%   10/27/2012        0.3             hme
%   04/16/2014        0.5             hme
%   07/18/2015        0.7             hme
%                

if nargin<2
  idx='';
end

N=length(v);
[~, n]=log2(N);
if (2^n-1)~=N
   error('Game has not the correct size!');
end

J=1:n;
S=bitset(0,J);
sS=S;
vi=v(S);
s_vi=vi*ones(n,1);
if s_vi>=v(N)
  error('Game is inessential!');
end
S(end+1)=N;
it=0:-1:1-n;
PlyMat=rem(floor(S(:)*pow2(it)),2);
%
% Defining imputation set.
%
A1=PlyMat(end,:);
A2=PlyMat(1:n,:);
B1=v(N);
B2=v(sS)';
% Defining the H polyhedron
H=struct('A',-[A1;A2],'B',-[B1;B2],'lin',(1:size(B1,1))');

% Calling cddmex
%V=cddmex('extreme',H);
% Trying to find the right position for plotting purpose.
if n<6
   P1 = Polyhedron( 'A', -A2, 'b', -B2, 'Ae', -A1, 'be', -B1);
   imp_vert=P1.V;
else
   V=cddmex('extreme',H);
   imp_vert=V.V;
end
if nargout == 2
   crst=[H.A,H.B];
elseif nargout > 2
  crst=[H.A,H.B];
  ip_v=imp_vert;
  if n>3
     if isempty(idx)
        y=range(ip_v);
        [~, idx]=min(y);
        ip_v(:,idx)=[];
     else
        ip_v(:,idx)=[];
     end
  else
     [X1,X2]=ToSimplex(ip_v);
     ip_v=[X1,X2];
  end
  P=Polyhedron(ip_v); % V-representation on R^n-1.
  ip_vol=volume(P);
end

