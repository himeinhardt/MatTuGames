function [lor_vert,lorst,lor_vol,P,Lsol]=LorenzSet(v,idx,tol)
% LORENZSET determines the Lorenz set of game v. 
% The cdd-library by Komei Fukuda is needed.
% It is recommended to install the cdd-library that accompanies
% the Multi-Parametric Toolbox 3.
% http://people.ee.ethz.ch/~mpt/3/
% 
%
% Usage: [lor_vert,lorst,lor_vol,P,Lsol]=LorenzSet(v,idx,tol) 
% Define variables:
%  output:
%  lor_vert     -- Matrix of the vertices of the Lorenz set. Output is numeric.
%  lorst        -- The core constraints.
%  lor_vol      -- Volume of the Lorenz set, if the Lorenz set is full dimensional,
%                  otherwise zero. 
%  P            -- Returns V- and H-representation (class Polyhedron)
%  Lsol         -- Returns structure variable of the Lorenz solution, see also the 
%                  function LorenzSol().
%  input:
%  v            -- A Tu-Game v of length 2^n-1. 
%  idx          -- Specifies the projection plane onto R^3. Default is the empty set.
%                  Will be computed internally.
%  tol          -- A positive tolerance value. Its default value is set to 10^9*eps.
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/12/2015        0.7             hme
%   01/14/2019        1.0             hme
%                

if nargin<2
  idx='';
  tol=10^6*eps;
elseif nargin < 3
  tol=10^8*eps;
end

N=length(v);
[~, n]=log2(N);
S=1:N;
PlyMat=zeros(N,n);
for k=1:n, PlyMat(:,k) = bitget(S,k);end
A1=PlyMat(1:N,:);
A1(N+1,:)=-PlyMat(end,:);

x1=ones(1,n)*v(N)/n;
crQ=CddCoreQ(v);
if crQ==0
   bQ=false;
else
   bQ=belongToCoreQ(v,x1);
end
Lsol=CPCore(v,x1,tol);
if bQ
   Lsol.Cp=x1;
   Lsol.D=0;
   Lsol.resid=zeros(1,n);
end
%% We can only get a single point of the Lorenz set by this appraoch!!!
adv=additive_game(Lsol.Cp);
adv=adv-tol;
adv(N)=v(N);

if crQ==1
   D1=-v(N);
   D2=adv';
% Defining the H polyhedron
   H=struct('A',A1,'B',[D2;D1],'lin',(1:size(D1,1))');

% Calling cddmex
   V=cddmex('extreme',H);
   lor_vert=V.V;
   if nargout == 2
      lorst=[H.A,H.B];
   elseif nargout > 2
      lorst=[H.A,H.B];
      cr_v=lor_vert;
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
      lor_vol=volume(P);
   end
else
     warning('Lor:Exit1','Core is empty! No Lorenz Set can be computed')
     lor_vert=[];
     lorst=[];
     lor_vol=[];
     P=[];
end
