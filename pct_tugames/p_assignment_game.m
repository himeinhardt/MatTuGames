function [v prof_mat]=p_assignment_game(sl_vec,prof_mat);
% P_ASSIGNMENT_GAME computes from an assignment problem (sl_vec,prof_mat) 
% the corresponding symmetric assignment game. If the problem is 
% not symmetric, it will be transformed into a symmetric one.
% This function needs a lot of memory for problems having 9 sellers onward.
% In this case at least 90.72 GB are needed. For 10 sellers at least
% 3628.8 GB are needed. Even for most modern 64-bit systems too much!!!
%
% Usage: [v prof_mat]=p_assignment_game(sl_vec,prof_mat)
% Define variables:
%  output:
%  v        -- An assignment game v of length 2^n-1.
%  prof_mat -- A square profit matrix.
%            
%  input:
%  sl_vec   -- A vector of sellers. From this vector, the vector of 
%              buyers will be constructed.
%              Example: sl_vec=[1 2 3 4 5];
%  prof_mat -- A non or square profit matrix.
%              use the function profit_matrix() to compute one.
%              or invoke
%              prfm=magic(5);
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   05/22/2011        0.1 alpha        hme
%   08/16/2012        0.2 beta         hme
%   03/02/2014        0.5              hme
%                

slnb=length(sl_vec);

if slnb>10
   error('Too much memory is requested! At least 3628.8 GB memory is needed to complete! Sorry!')
end

gr=max(size(prof_mat));
stm=size(prof_mat)==slnb;
if all(stm)==0
  warning('Assignment problem not symmetric. Transforming into a symmetric one.');
   sl_vec=1:gr;
   prof_mat(gr,gr)=0;
 else
end


m=length(sl_vec);
n=2*m;
by_vec=m+sl_vec;
N=2^n-1;
S=1:N;
v=zeros(1,N);


ps=sl_vec-1;
pb=by_vec-1;
sl=2.^ps;
by=2.^pb;
M1=sl*ones(m,1);
M2=by*ones(m,1);
sm1=SubSets(M1,n);
sm2=SubSets(M2,n);
Pm=cell(m);


sm=[sm1,sm2];
slC=setdiff(S,sm);
lC=length(slC);

parfor i=1:m
  for j=1:m
    Si=bitget(slC,i)==1;
    Sj=bitget(slC,by_vec(j))==1;
    Sij=Si & Sj;
    Pm{i,j}=prof_mat(i,j)*Sij;
  end
end


prm=perms(sl_vec);
spr=size(prm);
cb=spr(1);
A=speye(m);
Pk=false(m);
% This matrix needs for 8 sellers at least 2.52 GB memory.
v1=zeros(cb,lC); 


parfor k=1:cb
  Pk=A(prm(k,:),:);
  for i=1:m
   for j=1:m
     if Pk(i,j)==1
        v1(k,:)=Pm{i,j}+v1(k,:);
      else
     end
   end
  end
end

vm=max(v1);
v(slC)=vm;
