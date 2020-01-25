function PRepShap=p_replicate_Shapley(v,slc,tol,str)
% P_REPLICATE_SHAPLEY computes the game space that replicates the 
% Shapley value of a TU-game v from a linear basis using Matlab's PCT.
% For n>12 this function needs some time to complete.
%
% Usage: PRepShap=p_replicate_Shapley(v,scl,tol)
% Define variables:
% output:
%  shQ                   -- Indicates whether the games under consideration
%                           replicates the Shapley value. True (1)  or False (0).
%  shQv                  -- List of games that replicates the Shapley value.
%                           Vector of zeros/ones.                            
%  shm                   -- Solution matrix of computed Shapley values.
%  vm                    -- Game space spanned by the null space of the Shapley value.
%  str                   -- A string for getting a complete or reduced output.
%                           This useful to reduce the memory requirements for n>14. 
%                           Admissible values are:
%                           'red' that is, a reduced structure element is returned.
%                           'full' that is, the full output is returned (with matrices).
%                           Default is 'full' output.
%
%  input:
%  v      -- A Tu-Game v of length 2^n-1. 
%  scl    -- scaling factor, default is scl=1.
%  tol    -- tolerance value, default is tol=10^8*eps.
%

%
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   06/22/2013        0.4             hme
%

if nargin < 1
   error('At least the game must be defined!')
elseif nargin < 2
  slc=1;
  tol=10^8*eps;
  str='full';
elseif nargin < 3
  tol=10^8*eps;
  str='full';
elseif nargin < 4
  str='full';
end


N=length(v);
[~, n]=log2(N);


slb=p_linear_basis(n,'sparse');
% To release memory for n>15
k=1:n;
sC=2.^(k-1);
slb(sC,:)=[];
[d1, ~]=size(slb);
vm=repmat(v,d1,1)+slc*slb;
clear slb;
ed=d1+1;
vm(ed,:)=v;
shm=zeros(ed,n);
parfor k=1:ed
    shm(k,:)=ShapleyValue(vm(k,:));
end
shmE=repmat(shm(ed,:),ed,1)';
shQv=all(abs(shm'-shmE)<tol);
shQ=all(shQv);
if strcmp(str,'full')
   PRepShap=struct('shQ',shQ,'shQv',shQv,'shm',shm,'vm',vm);
elseif strcmp(str,'red')
   PRepShap=struct('shQ',shQ,'shQv',shQv);
else
   PRepShap=struct('shQ',shQ,'shQv',shQv);
end
