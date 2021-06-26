function [spU,ocU]=DecomposeVN(n,method)
% DECOMPOSEVN computes a subspace and its orthogonal complement of the game space VN (G^n) for n-players.
% Requires spnull to compute the null space of a sparse matrix.
% Source: Beal et al. (2013) Proposition 4.
%
% Usage: [spU,ocU]=DecomposeVN(n)
%
% Define variables:
%  output:
%  spU      -- Cell output of size (n,1). Cell element i contains the
%              spanning system of the subspace U of player i. 
%  ocU      -- Cell output of size (n,1). Cell element i contains
%              the spanning system of the orthogonal complement of spU.
%
%  input:
%  n        -- An integer greater or equal to 2.
%  method   -- A string to format the matrices. Permissible methods
%              to format the matrices are 'full','sparse' or the
%              empty string '', to invoke the default, which is sparse. 
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
%   02/18/2014        0.5             hme
%
if nargin < 2
   method='sparse';
else
   if isempty(method)
      method='sparse'; 
   end 
end    

N=2^n-1;
S=1:N;
if strcmp(method,'sparse')
   dr=speye(N);
else
   dr=eye(N); 
end    
int=0:-1:1-n;
cmat=(rem(floor(S(:)*pow2(int)),2)==1)';
np=sum(cmat,1);
spU=cell(n,1);
ocU=cell(n,1);
for pl=1:n
    npl=S(bitget(S,pl)==0);
    mnpl=dr(npl,:);
    for sz=1:n
       sm=S(np==sz);
       pm=sm(bitget(sm,pl)==1);
       vl=length(pm);
       msz(sz,:)=ones(1,vl)*dr(pm,:);
    end
    if strcmp(method,'sparse') 
       spU{pl}=[mnpl;msz];
       try
         ocU{pl}=spnull(spU{pl});
       catch
         ocU{pl}=null(spU{pl});
       end
    else
       spU{pl}=[mnpl;msz];
       ocU{pl}=null(spU{pl});
    end    
end
