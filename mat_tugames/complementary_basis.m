function cb=complementary_basis(n,method)
% COMPLEMENTARY_BASIS computes a basis of the n-person TU game space.
% 
% Usage: gb=complementary_basis(n,method)
%

% Define variables:
%  output:
%  cb       -- The complementary basis of an n-person TU game space. 
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
%   02/16/2022        1.9.1           hme
%                

if nargin < 2
   method='full';
else
   if isempty(method)
      method='full';
   end
end

N=2^n-1;
N1=N-1;
cb=false(N);
S=1:N;
int=1-n:1:0;

for k=1:N1
   aS=bitand(S,k)==0;
   cb(k,:)=aS;
end
cb=[true(1,N); cb]; %% as(N) is as(0)=\vec{1}.
cb(N+1,:)=[];
if strcmp(method,'sparse')
   cb=sparse(cb);
end

