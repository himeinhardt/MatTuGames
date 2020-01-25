function [F,dfm]=MLextension(clv)
% MLEXTENSION computes the multi-linear extension of a TU-game v 
% using the multi-linear extension. It requires Matlab's Symbolic 
% Math Toolbox.   
%
% Usage: F=clv.MLextension
%
% Define variables:
%  output:
%  F        -- The multi-linear extension of game v.
%  dfm      -- The matrix of partial derivatives.
%              Must be integrated from 0 up 1 to get
%              the Shapley value.
%
%  input:
%  clv      -- TuGame class object.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   09/29/2014        0.5             hme
%

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;

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
  dfm{jj}=dfy;
end

