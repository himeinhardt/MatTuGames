function [x1, fmin]=cplex_cfr_nucl(clv,F,tol)
% CPLEX_CFR_NUCL computes the nucleolus of game v with coalition fromation restrictions 
% using cplexmex.
%
% Source: Granot et al. (1978), Characterization sets for the nucleolus. IJGT.
%
% http://www-01.ibm.com/software/websphere/ilog/
% (compatible with CPLEX Version 12.10.0 and higher)
% 
%
% Usage: [x, fmin]=clv.cplex_cfr_nucl(F)
% Define variables:
%  output:
%  x1        -- The nucleolus of game vF.
%  fmin      -- The minmax excess value.
%
%  input:
%  clv      -- TuGame class object.
%  F        -- For instance, a characterization set for the nucleolus.
%              F must contain the grand coalition N.
%  tol      -- Tolerance value. Its default value is set to 10^8*eps.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/15/2017        0.9             hme
%   04/04/2020        1.9             hme
%                


warning('off','all');
if nargin<2
   error('A collection of sets F is required!');
elseif nargin<3
 tol=10^8*eps;
end


v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
essQ=clv.tuessQ;
si=clv.tuvi;
if essQ==0
   error('Game is not essential!')
end

%% Addining the singleton coalitions to F, 
%% since v is definded over F ans si.
F=unique([F,si]);
lf=length(F);

% solver parameter
ra = clv.reasonable_outcome();
cvr=vi==ra;
if any(cvr)
   fi=find(cvr);
   ra(fi)=Inf;
end
lb=[vi,-Inf]';
ub=[ra,Inf]';

x0=[];
mtv=verLessThan('matlab','9.1.0');
if mtv==1
  options = cplexoptimset('MaxIter',128,'Simplex','on','Display','off');
else 
%  options = cplexoptimset('MaxIter',128,'Algorithm','primal','Display','off');
  options.largescale='on';
  options.algorithm='dual-simplex';
  options.tolfun=1e-10;
  options.tolx=1e-10;
  options.tolrlpfun=1e-10;
  %%%% for dual-simplex
  % opts.MaxTime=9000;
  options.preprocess='none';
  options.tolcon=1e-6;
  options.maxiter=128;
  options.display='off';
  options.threads=3;
end

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
for k=1:n, A1(:,k) = -bitget(F,k);end
A1(end+1,:)=-A1(end,:);
A1(:,end+1)=-1;
A1(lf:lf+1,end)=0;
A2=sparse(A1);
B1=[-vF';vF(lf)];
C=[zeros(n,1);1];

while 1
  [xmin,fmin,exitflag,~,lambda]=cplexlp(C,A2,B1,[],[],lb,ub,x0,options);
  if exitflag ~= 1     
     warning('on','all');
     warning('Prn:Exit','Probably no nucleolus found!')
     break;
  end
  x=xmin;
  x1=x';
  x1(end)=[];
  bS1=(find(lambda.ineqlin>tol))';
  bS1(end)=[];
  bA=find(A1(:,end)==0)';
  bS2=setdiff(bS1,bA);
  if isempty(bS2)==1
     warning('on','all');
     break;
  end
  bA(end)=[];
  bA=F(bA);
  it=0:-1:1-n;
  mS2=rem(floor(bA(:)*pow2(it)),2);
  tmS2=mS2';
  rk=rank(mS2);
  ov=ones(1,n);
  wgh=pinv(tmS2)*ov';
  posQ=all(wgh>-tol);
  if rk==n && posQ == 1
     warning('on','all');
     break;
  end
  A1(bS2,end)=0;
  A2=sparse(A1);
  B1(bS2)=B1(bS2)+fmin;
end
