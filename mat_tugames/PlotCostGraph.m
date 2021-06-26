function PlotCostGraph(cm,mE)
% PLOTCOSTGRAPH plots from a cost matrix the associated cost spanning graph. 
%
% Usage: PlotCostGraph(cm)
%
% Define variables:
%  output:
%           -- A plot of the cost spanning graph.
%            
%  input:
%
%  cm       -- A square cost matrix (n+1xn+1) derived from a cost
%              spanning graph problem. For instance, for a four
%              person game the size of the matrix must be (5x5). The source
%              is player 1.
%  mE       -- An alternative labeled edge matrix like 
%              mE=[0 2 4;4 4 5] instead of [0 1 2;2 2 3] got by 
%              internal operations.
%

%
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   04/21/2021        1.9             hme
%

if nargin < 2
   mE=[];
elseif nargin == 2
   s2=mE(1,:);
   t2=mE(2,:);
end	

[m1,m2]=size(cm);
if m1~=m2
   error('Matrix is not square');
else
   n=m1;
end
N=2^n-1;
n1=n-1;
N1=2^n1-1;

upe=true(n);
ann=n^2;
ar=1:ann;
szA=[n,n];
A=reshape(ar,szA);
UA=triu(upe,1);
ind=A(UA)';
wghs=cm(UA)';
slc=wghs>0;
wghs=wghs(slc);
[s,t]=ind2sub(szA,ind);
s=s(slc);
t=t(slc);
if isempty(mE)
   s1=s-1;
   t1=t-1;
else
   s1=s2;
   t1=t2;
end
s0=string(s1);
t0=string(t1);
G = graph(s0,t0,wghs);
%%
%% Plotting the cost graph
clf;
figure(1); hold on
gr = plot(G,'EdgeLabel',G.Edges.Weight);
box on;
set(gca,'XTick',[],'YTick',[]);
%box off;
hold off
