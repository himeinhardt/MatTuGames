function v=p_PermutationGame(asm)
% P_PERMUTATIONGAME computes from an assignment matrix the corresponding 
% permutation game using Matlab's PCT. 
% Requires for n>10 more than 64 GB physical memory.
%
% Usage: v=p_PermutationGame(asm)
%
% Define variables:
%  output:
%  v        -- A permutation game v of length 2^n-1.
%            
%  input:
%
%  asm      -- A square assignment matrix.
%              To construct one invoke, for instance,
%              asm=magic(5);
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   09/20/2014        0.5             hme
%                

[n1,n2]=size(asm);
if n1~=n2
   error('Matrix is not square!'); 
end 
n=n1;
pl_vec=1:n;
N=2^n-1;
S=1:N;
int=0:-1:1-n;
mat=(rem(floor(S(:)*pow2(int)),2)==1);
v=zeros(1,N);
A=eye(n);

parfor ss=1:N
 pS=pl_vec(mat(ss,:));
 prm=perms(pS);
 [spm1,spm2]=size(prm);
 pmat=A(pS,:)==1;
 ak=zeros(1,spm1);
 for k=1:spm1
    pt=prm(k,:); 
    pfrd=asm(pt,:);
    ak(k)=sum(pfrd(pmat));
 end
 v(ss)=max(ak);
end
    
