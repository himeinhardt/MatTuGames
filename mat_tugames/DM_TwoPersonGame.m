function [v_ij,A] = DM_TwoPersonGame(v,x)
% DM_TWOPERSONGAME computes from (v,x) all reduced two-person games at x of
% game v.
% 
% Usage: [v_ij,A] = DM_TwoPersonGame(v,x)
%
% Define variables:
%  output:
%  v_ij     -- All Davis-Maschler reduced two-person games w.r.t. x.
%  A        -- Matrix of standard solutions.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n). Must be efficient.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   10/10/2015        0.7             hme
%

N=length(v);
[~, n]=log2(N);
S=1:N;
k=1:n;
sC=2.^(k-1);
A=zeros(n);
lt=(n+1)*n/2;
v_ij=zeros(lt,3);
ct=1;
ad_x=additive_game(x);
for ii=1:n-1
  for jj=ii+1:n
      Q=S(bitget(S(bitget(S,ii)==0),jj)==0);
      sS2=bitor(sC(ii),Q);
      sS3=bitor(sC(jj),Q);
      A(ii,jj)=max([v(sS2)-ad_x(Q),v(sC(ii))]);
      A(jj,ii)=max([v(sS3)-ad_x(Q),v(sC(jj))]);
      vN=A(ii,jj)+A(jj,ii);
      v_ij(ct,:)=[A(ii,jj),A(jj,ii),vN];
      ct=ct+1;
  end 
end 
