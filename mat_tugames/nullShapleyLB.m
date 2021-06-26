function [nSh,CnSh]=nullShapleyLB(n,method)
% NULLSHAPLEYLB determines a basis of the null space for the
% Shapley-value for n-persons by a linear basis approach.
% Source: Yokote et al. (2013)
%
% Usage: kSh=nullShapleyLB(n,'sparse')
%
% Define variables:
%  output:
%  nsh      -- A basis of the null space for the Shapley value
%              for n-persons.    
%  CnSh     -- A basis of the complement space.
%
%  input:
%  n        -- Number of players involved.
%  method   -- A string to format the matrices. Permissible methods
%              to format the matrices are 'full','sparse' or the
%              empty string '', to invoke the default, which is sparse.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/28/2014        0.5             hme
%                
    

if nargin < 1
   error('At least the number of players must be defined!')
elseif nargin < 2
   method='full';
end



slb=linear_basis(n,method);
k=1:n;
sC=2.^(k-1);
CnSh=slb(sC,:);
slb(sC,:)=[];
nSh=slb;
