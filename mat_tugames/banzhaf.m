function [bzv theta]=banzhaf(sv);
% BANZHAF computes the Banzhaf/Coleman index of a simple game sv.
%
% Usage: [bzv theta]=banzhaf(sv)
% Define variables:
%  output:
%  bzv       -- The Banzhaf/Coleman index (normalized).
%  theta     -- A vector of swings.
%
%  input:
%  sv        -- A simple game of length 2^n-1.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/15/2010        0.1 beta        hme
%   06/27/2012        0.2 beta        hme
%   10/27/2012        0.3             hme
%   03/08/2014        0.5             hme
%                



N=length(sv);
[~, n]=log2(N);
S=1:N;
bzv=zeros(1,n);
theta=zeros(1,n);

for k=1:n
  plk=bitget(S,k)==1;
  nplk=plk==0;
  Clk=S(plk);
  Clnk=S(nplk);
  w_clk=sv(Clk);
  nw_clk=sv(Clnk);
  sum_wk=w_clk*ones(length(w_clk),1);
  sum_nwk=nw_clk*ones(length(nw_clk),1);
  theta(k)=sum_wk-sum_nwk;
end


bzv=theta/(theta*ones(n,1));
