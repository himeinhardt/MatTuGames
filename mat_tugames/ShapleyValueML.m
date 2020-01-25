function [shv,F,dfm]=ShapleyValueML(v,str)
% SHAPLEYVALUEML computes the Shapley-value of a TU-game v 
% using the multi-linear extension. It requires Matlab's Symbolic 
% Math Toolbox.   
%
% Credits: We are thankful to Amro from StackOverFlow correcting
%          a bug.
%
% Usage: [shv,F,dfm]=ShapleyValueML(v,'sym')
%
% Define variables:
%  output:
%  shr      -- The Shapley-value of a TU-game v (rational numbers). 
%  F        -- The multi-linear extension of game v.
%  dfm      -- The matrix of partial derivatives.
%              Must be integrated from 0 up 1 to get
%              the Shapley value.
%
%  input:
%  v        -- A TU-Game of length 2^n-1.
%  str      -- Converts symbolic expression into numeric form.
%              Permissible strings are:
%              'num' numeric form.
%              'sym' symbolic form.
%              ''    empty string, which is symblic form (default).
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   09/27/2014        0.5             hme
%

if nargin <2
   str='sym';
end

N=length(v);
[~, n]=log2(N);
if (2^n-1)~=N
    error('Game has not the correct size!');
end
if N==1
  shv=v;return;
elseif iscolumn(v)
    v=v';
end

S=1:N;
it=0:-1:1-n;
mat=(rem(floor(S(:)*pow2(it)),2)==1);
cmat=mat==0;
x=sym('x',[1 n]);
mx=1-x;
y = sym('y');
vy=ones(1,n)*y;
F=0;
if strcmp(str,'num')
   shv=zeros(1,n);
else
   shv=sym('shv',[1,n]);
end
dfm=cell(1,n);

for ss=1:N
  pd1=x(mat(ss,:));
  pd2=mx(cmat(ss,:));
  pd=prod(pd1)*prod(pd2)*v(ss);
  F=pd+F;
end
F=expand(F);

for jj=1:n
  dF=diff(F,x(jj));
  dfy=subs(dF,x,vy);
  if strcmp(str,'num')
     mf=matlabFunction(dfy);
     shv(jj)=integral(mf,0,1);
  else
     shv(jj)=int(dfy,0,1);
  end
  dfm{jj}=dfy;
end

