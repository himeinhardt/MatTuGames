function [bzv theta]=p_banzhaf(sv)
% P_BANZHAF computes the Banzhaf/Coleman index of a simple game sv.
% Using Matlab's PCT.
%
% Usage: [bzv theta]=p_banzhaf(sv)
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
%   05/21/2011        0.1 alpha        hme
%   06/27/2012        0.2 beta         hme
%   10/27/2012        0.3              hme
%   05/16/2014        0.5              hme
%                



N=length(sv);
[~, n]=log2(N);
S=1:N;
bzv=zeros(1,n);
theta=zeros(1,n);

parfor k=1:n
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
