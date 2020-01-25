function lb=linear_basis(n,method)
% LINEAR_BASIS computes the linear basis of the n-person TU game space.
% 
% Usage: lb=linear_basis(n,method)
%

% Define variables:
%  output:
%  lb       -- The linear basis of an n-person TU game space. 
%
%  input:
%  n        -- number of players involved.
%  method   -- A string to format the matrix. Permissible methods
%              to format the matrices are 'full','sparse' or the
%              empty string '', to invoke the default, which is full.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   06/21/2013        0.4             hme
%   02/15/2014        0.5             hme
%                

if nargin < 2
   method='full';
else
   if isempty(method)
      method='full';
   end
end

N=2^n-1;
lb=false(N);
S=1:N;
int=1-n:1:0;

for k=1:N
   aS=bitand(S,k);
   mat=rem(floor(aS(:)*pow2(int)),2)';
   lb(k,:)=ones(1,n)*mat==1;
end

if strcmp(method,'sparse')
   lb=sparse(lb);
end
