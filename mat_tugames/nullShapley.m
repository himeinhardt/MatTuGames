function [kSh,CkSh]=nullShapley(n)
% NULLSHAPLEY determines a basis of the null space for the
% Shapley-value for n-persons.
% Source: Faigle and Grabisch (2014)
%
% Usage: kSh=nullShapley(n)
%
% Define variables:
%  output:
%  ksh      -- A basis of the null space for the Shapley value
%              for n-persons.    
%  CkSh     -- A basis of the complement space.
%
%  input:
%  n        -- Number of players involved.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   03/08/2014        0.5             hme
%                
    
gb=game_basis(n);
N=2^n-1;
k=1:n;
pl=2.^(k-1);
CkSh=gb(:,pl);
int=1-n:1:0;
S=1:N;
S(pl)=[];
mat=rem(floor(S(:)*pow2(int)),2)==1; 
cS=mat*ones(n,1);
clear mat;
lS=N-n;
kSh=zeros(lS,N);
for ii=1:lS
  T=S(ii);
  ci=k(bitget(T,k)==1);
  lci=length(ci);
  kSh(ii,:)=gb(:,T)-(1/cS(ii))*CkSh(:,ci)*ones(lci,1);
end
