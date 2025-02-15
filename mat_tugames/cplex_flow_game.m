function [v,x]=cplex_flow_game(E,c,ow)
% CPLEX_FLOW_GAME computes from a flow problem (E,c,ow) a TU flow game using cplexmex.
%
% Usage: v=cplex_flow_game(E,c,ow)
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
%  v=cplex_flow_game(E,c,ow)
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
%   02/24/2018        0.9             hme
%   04/04/2020        1.9             hme
%
warning('off','all');
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
% lower boundary
lb=zeros(m,1);
ub=c';
%
% Options setting
% solver parameter
mtv=verLessThan('matlab','9.1.0');
if mtv==1
  opts = cplexoptimset('MaxIter',128,'Simplex','on','Display','off');
else
%  opts = cplexoptimset('MaxIter',128,'Algorithm','primal','Display','off');
  opts.largescale='on';
  opts.algorithm='dual-simplex';
  opts.tolfun=1e-10;
  opts.tolx=1e-10;
  opts.tolrlpfun=1e-10;
  %%%% for dual-simplex
  % opts.MaxTime=9000;
  opts.preprocess='none';
  opts.tolcon=1e-6;
  opts.maxiter=128;
  opts.display='off';
  opts.threads=3;
end
%
% Matrix construction
A2=eye(m);
B1=zeros(b,1);
B2=c';
[xmax, fmax, status, extra] = cplexlp(-C,A2,B2,A1,B1,lb,ub,[],opts);
v(N)=C'*xmax;
x(N,:)=xmax';


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
   B3(idx2)=0;
   % objective function for coalition S.
   l1=size(Es2);
   ls=1:l1;
   zc2=ls(Es2(:,2)==b+1)';
   zc3=zc(ismember(E(zc,1),Es2(zc2,1)));
   C=zeros(m,1);
   C(zc3)=1;
   [xmax, fmax, status, extra] = cplexlp(-C,A2,B3,A1,B1,lb,ub,[],opts);
   v(S)=C'*xmax;
   x(S,:)=xmax';
end
warning('on','all');
