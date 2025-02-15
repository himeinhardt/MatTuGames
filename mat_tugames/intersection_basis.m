function gb=intersection_basis(n,method)
% INTERSECTION_BASIS computes a basis of the n-person TU game space.
% 
% Usage: gb=intersection_basis(n,method)
%

% Define variables:
%  output:
%  gb       -- The intersection basis of an n-person TU game space. 
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
%   06/24/2021        1.9             hme
%                

if nargin < 2
   method='full';
else
   if isempty(method)
      method='full';
   end
end

N=2^n-1;
gb=false(N);
S=1:N;
int=1-n:1:0;

for k=1:N
   aS=bitand(S,k)>0;
   gb(k,:)=aS;
end

if strcmp(method,'sparse')
   gb=sparse(gb);
end
