function lb=PS_GameBasis(n,method)
% PS_GAMEBASIS computes the basis for the class of PS games.
% 
% Usage: lb=PS_GameBasis(n,method)
%
% Resource: Funaki and Yokote (2020), Chap 7, Handbook of the Shapley value,  eds. E. Algaba, V. Fragnelli and J. Sanchez-Soriano.
%
% Define variables:
%  output:
%  lb       -- The basis of the PS game space. 
%
%  input:
%  n        -- Number of players involved.
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
%   04/15/2020        1.9             hme
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
pl=1:n;
dx=false(1,N);
gb=game_basis(n);

for k=1:N
   aS=bitand(S,k);
   mat=rem(floor(aS(:)*pow2(int)),2)';
   T=nnz(bitget(k,pl));
   if T==1;
     lb(:,k)=gb(:,k);
   elseif ~mod(T,2)
     lb(:,k)=ones(1,n)*mat==T/2;
   else 
     dx(k)=true;	   
   end
end
lb(:,S(dx))=[];

if strcmp(method,'sparse')
   lb=sparse(lb);
end

