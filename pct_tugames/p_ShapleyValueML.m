function [shv,F,dfm]=p_ShapleyValueML(v)
% SHAPLEY_VALUE computes the Shapley-value of a TU-game v 
% using the multi-linear extension. It requires Matlab's PCT and 
% and the Symbolic Math Toolbox. But the Symbolic Math Toolbox
% does not work together with the PCT. 
%
% Usage: [shv,F,dfm]=p_ShapleyValueML(v)
%
% Define variables:
%  output:
%  shv      -- The Shapley-value of a TU-game v.
%  F        -- The multi-linear extension of game v.
%  dfm      -- The matrix of partial derivatives.
%              Must be integrated from 0 up 1 to get
%              the Shapley value.
%
%  input:
%  v        -- A TU-Game of length 2^n-1.
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


N=length(v);
[~, n]=log2(N);
S=1:N;
int=0:-1:1-n;
mat=(rem(floor(S(:)*pow2(int)),2)==1);
cmat=mat==0;
x=sym('x',[1 n]);
mx=1-x;
y = sym('y');
vy=ones(1,n)*y;
F=0;
shv=zeros(1,n);
dfm=cell(1,n);

parfor ss=1:N
  pd1=x(mat(ss,:));
  pd2=mx(cmat(ss,:));
  pd=prod(pd1)*prod(pd2)*v(ss);
  F=pd+F;
end
F=expand(F);

parfor jj=1:n
  dF=diff(F,x(jj));
  dfy=subs(dF,x,vy);
%% Does not work!! Matlab bug??? 
%  mf=matlabFunction(dfy);
%  shv(jj)=integral(mf,0,1);
%%
%% The best would be to use:
%%
%   shv(jj)=int(dfy,0,1)
%% but it cannot be used in a function.
  dfm{jj}=dfy;
end

