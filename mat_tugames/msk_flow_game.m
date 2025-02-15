function [v,A,B,C,x]=msk_flow_game(E,c,ow)
% MSK_FLOW_GAME computes from a flow problem (E,c,ow) a TU flow game using mosekmex.
% 
% MSK-SOLVER: http://www.mosek.com/
%
% Usage: v=msk_flow_game(E,c,ow)
%
% Define variables:
% output:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- A matrix of optimal solution vectors for all
%              coalition S.
% input: 
%  E        -- An edge matrix of size (lx2) or a cell of numel l.
%              The source must be given by 0, and the sink by the
%              number of vertices plus one. However, the edge matrix
%              can also be of size (lx3) then c can be empty. This is
%              due to be compatible with the GrTheory toolbox. In this
%              case there is no need to adjust the edge matrix from GrTheory.
%  c        -- A capacity vector of max flow for each edge. Can be set
%              to [] if the edge matrix has size(lx3). The capacity
%              vector will then be obtained from E(:,3).
%  ow       -- An ownership vector of the edges.
%
%
% Example:
% Define a matrix of edges given by
% E =
%   0   0   1   1   1   4   2   3   3   4   5
%   1   2   2   3   4   2   5   4   6   5   6
%
% or equivalently
%
% E={[0 1];[0 2];[1 2];[1 3];[1 4];[4 2];[2 5];[3,4];[3,6];[4,5];[5,6]};
% Do not forget to set the separator by semicolon (;) not by comma.
%
% We have 11 edges here. Furthermore, in total we have 5 vertices, and one
% source given by number 0, and a sink by number 6.
% Then specify the capacity vector of max flow of each edge.
% c =
%   2   2   1   2   1   1   2   1   2   1   3
%
% Finally, determine to which player the edges are owned
% ow =
%   1   2   3   1   1   2   2   2   3   3   1
%
% Here, we have a player set of {1,2,3}, hence n=3.
%
% Now, invoke
%  v=msk_flow_game(E,c,ow)
% to get
% v =
%     0     0     2     0     2     0     4
%
%


%
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   02/10/2016        0.8             hme
%

if iscell(E)
   E=cell2mat(E);
else 
   [m,cs]=size(E);
   if m==2
      E=E';
   elseif m==3 % to be compatible with the GrTheory toolbox.
      c=E(3,:);
      E(3,:)=[];
      E=E';
   end
   if cs==3  % to be compatible with the GrTheory toolbox.
      c=E(:,3)';
      E(:,3)=[];
   end   
end
[m,~]=size(E);
n=max(ow);
N=2^n-1;
v=zeros(1,N);
x=zeros(N,m);
l=1:m;
% to be compatible with the GrTheory toolbox.
% if the source is given by a positive integer.
me=min(E);
if me(1)>0
   E=E-me(1);
end
vn=max(E);
ig=cell(1,n);
og=cell(1,n);

% Max flow for the grand coalition.
b=vn(1);
A1=zeros(b,m);
for ii=1:b
    og=l(E(:,1)==ii);
    ig=l(E(:,2)==ii);
    A1(ii,og)=-1;
    A1(ii,ig)=1;    
end
zc=l(E(:,2)==b+1)';
% objective function
C=zeros(m,1);
C(zc)=1;
%
% Changing parameter values to increase precision.
[rcode,res] = mosekopt('param echo(0)');
param=res.param;
%param.MSK_IPAR_INTPNT_BASIS   = sc.MSK_OFF;
%param.MSK_DPAR_INTPNT_TOL_REL_GAP = 1.0000e-12; % Adjust this value if the solution is not correct.
%param.MSK_IPAR_OPTIMIZER = 5;  % Using dual simplex. MSK 7
param.MSK_IPAR_OPTIMIZER ='MSK_OPTIMIZER_DUAL_SIMPLEX'; % MSK 8
%param.MSK_DPAR_BASIS_TOL_X = 1.0e-9;
%param.MSK_DPAR_BASIS_TOL_S = 1.0e-9;

prob.c=C;
% lower/upper boundary
prob.blx=zeros(m,1);
prob.bux=c';

%
% Matrix construction
A2=eye(m);
B1=zeros(b,1);
B2=c';
A=[A1;A2];
Bu=[B1;B2];
B2l=-inf(m,1);
Bl=[B1;B2l];
prob.a=sparse(A);
prob.buc=Bu;
prob.blc=Bl;

[rcode,res] = mosekopt('maximize echo(0)',prob,param);
sol=res.sol;
xx=sol.bas.xx;
v(N)=C'*xx;
x(N,:)=xx';


% Flows for sub-coalitions.
it=0:-1:1-n;
k=1:n;
for S=1:N-1
   ci=rem(floor(S(:)*pow2(it)),2)==1;
   rS=ismember(ow,k(ci));
   idx=l(rS);
   idx2=l(ismember(l,idx)==0);
%   Es=E(idx2,:);
   Es2=E(idx,:);
   % Matrix and boundary vector for coalition S.
   B3=B2;
   B3l=B2l;
   B3(idx2)=0;
   B3l(idx2)=0;
   B3u=[B1;B3];
   B3l=[B1;B3l];
   prob.buc=B3u;
   prob.blc=B3l;
   c2=c';
   c2(idx2)=0;
   prob.bux=c2;
   % objective function for coalition S.
   l1=size(Es2);
   ls=1:l1;
   zc2=ls(Es2(:,2)==b+1)';
   zc3=zc(ismember(E(zc,1),Es2(zc2,1)));
   C=zeros(m,1);
   C(zc3)=1;
   prob.c=C;
   [rcode,res] = mosekopt('maximize echo(0)',prob,param);
   sol=res.sol;
   xx=sol.bas.xx;
   v(S)=C'*xx;
   x(S,:)=xx';
end
