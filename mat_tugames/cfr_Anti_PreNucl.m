function [x1, fmin]=cfr_Anti_PreNucl(v,F,tol)
% CFR_ANTI_PRENUCL computes the anti pre-nucleolus of game v with coalition formation restrictions
% using the optimization toolbox.
% Uses now Dual-Simplex (Matlab R2015a).
%
% Usage: [x, fmin]=cfr_Anti_PreNucl(v,F,tol)
% Define variables:
%  output:
%  x1        -- The anti pre-nucleolus of game v w.r.t. F.
%  fmin      -- The maxmin excess value.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1.
%  F        -- For instance, a characterization set for the nucleolus.
%              F must contain the grand coalition N. 
%  tol      -- Tolerance value. Its default value is set to 10^8*eps.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/09/2017        0.9             hme
%   05/01/2024        1.9.2           hme
%                


if nargin<2
   error('A collection of sets F is required!');
elseif nargin<3
 tol=10^6*eps;
end


N=length(v);
[~, n]=log2(N);
k=1:n;
si=bitset(0,k);
%% Addining the singleton coalitions to F, 
%% since vF is definded over F and si.
F=unique([F,si]);
lf=length(F);


%% F should contain the grand coalition for defining vF.
S=1:N;
lfNq=F(end)~=N;
if lfNq
   F(end+1)=N;
   lf=lf+1;
end
CS=S(ismember(S,F)==0);
vF=v;
vF(CS)=[];

if lfNq
   vF(end)=0;
end



for k=1:n, A1(:,k) = bitget(F,k);end
A1(lf+1,:)=-A1(end,:);
A1(:,end+1)=1;
A1(lf:lf+1,end)=0;
A2=sparse(A1);
B1=[vF';-vF(lf)];
C=[zeros(1,n),-1];

%ra = smallest_amount(v);
%ub=[-ra,inf];
x1=[];
lb=[-inf(1,n),-inf];
%ub=[ra,Inf];
ub=[];
opts.Simplex='on';
opts.Display='off';
%opts.ActiveSet='on';
opts.LargeScale='on';
mth1=verLessThan('matlab','24.1.0');
if mth1==0,
    opts.Algorithm='dual-simplex-highs';
else
    opts.Algorithm='dual-simplex';
end
opts.TolFun=1e-10;
opts.TolX=1e-10;
opts.TolRLPFun=1e-10;
%opts
%% for dual-simplex
opts.MaxTime=9000;
opts.Preprocess='none';
opts.TolCon=1e-6;
opts.MaxIter=10*(N+n);

while 1
  try
    [xmin,fmin,exitflag,~,lambda]=linprog(C,A2,B1,[],[],lb,ub,opts);
  catch %% old api (before R2022a) with initial value.
   [xmin,fmin,exitflag,~,lambda]=linprog(C,A2,B1,[],[],lb,ub,[],opts);
  end	
  x=xmin';
  x1=x;
  if isempty(x1) == 1
     warning('Prn:Exit0','Probably no pre-nucleolus found!')
     break;
  end
  x1(end)=[];
  bS1=find(lambda.ineqlin'>tol);
  bA=find(A1(:,end)==0)';
  bS2=setdiff(bS1,bA);
  if isempty(bS2)==1
     break;
  end
  it=0:-1:1-n;
  bA(end)=[];
  bA=F(bA);
  mS2=rem(floor(bA(:)*pow2(it)),2);
  tmS2=mS2';
  rk=rank(mS2);
  ov=ones(1,n);
  wgh=pinv(tmS2)*ov';
  posQ=all(wgh>-tol);
  if exitflag ~= 1
     warning('Prn:Exit1','Probably no pre-nucleolus found!')
     break; 
  elseif rk==n && posQ == 1
     x1=full(x1); 
     break;
  end
  A1(bS2,end)=0;
  A2=sparse(A1);
  B1(bS2)=B1(bS2)+fmin;
  y=x1;
end
