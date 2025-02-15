function [bpv bzv]=BanzhafPenrose(sv);
% BANZHAFPENROSE computes the Banzhaf/Penrose and Banzhaf/Coleman index of a simple game sv.
%
%
% Usage: [bpv bzv]=BanzhafPenrose(sv)
%
% Define variables:
%  output:
%  bpv       -- The Banzhaf/Penrose index (normalized Banzhaf/Penrose index to one).
%  bzv       -- The Banzhaf/Coleman index (normalized Banzhaf/Coleman index by N).
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
%   05/25/2014        0.5             hme
%   09/09/2021        1.9.1           hme
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
bzv=theta/N; % Coleman
bpv=bzv./(theta*ones(n,1)); % Penrose
