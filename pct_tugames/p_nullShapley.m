function [kSh,CkSh]=p_nullShapley(n)
% P_NULLSHAPLEY determines a basis of the null space for the
% Shapley-value for n-persons using Matlab's PCT.
% Source: Faigle and Grabisch (2014)
%
% Usage: [kSh,CkSh]=p_nullShapley(n)
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
    
gb=p_game_basis(n);
N=2^n-1;
k=1:n;
pl=2.^(k-1);
CkSh=gb(:,pl);
int=1-n:1:0;
S=1:N;
parfor jj=1:n
   mat(:,jj)=bitget(S,jj)==1;
end
mat(pl,:)=[];
S(pl)=[];
cS=mat*ones(n,1);
clear mat;
lS=N-n;
kSh=zeros(lS,N);
parfor ii=1:lS
  T=S(ii);
  ci=k(bitget(T,k)==1);
  lci=length(ci);
  kSh(ii,:)=gb(:,T)-(1/cS(ii))*CkSh(:,ci)*ones(lci,1);
end
